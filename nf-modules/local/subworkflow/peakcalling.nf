/***********************
 * Peak calling
 */

include { samtoolsFlagstat } from '../../common/process/samtoolsFlagstat'
include { macs2 as macs2Sharp } from '../process/macs2' 
include { macs2 as macs2Broad } from '../process/macs2'
include { epic2 } from '../process/epic2'
include { frip } from '../process/frip'
include { annotPeaksHomer } from '../process/annotPeaksHomer'
include { peakQC } from '../process/peakQC'
include { IDR } from '../process/IDR'

Channel
  .fromPath("$baseDir/assets/peak_count_header.txt")
  .set { chPeakCountHeader }
Channel
  .fromPath("$baseDir/assets/frip_score_header.txt")
  .set { chFripScoreHeader }
Channel
  .fromPath("$baseDir/assets/peak_annotation_header.txt")
  .set{ chPeakAnnotationHeader }

// Create special channel to deal with no input cases
Channel
  .from( ["NO_INPUT", file("NO_FILE_BAM"), file("NO_FILE_BAI")] )
  .toList()
  .set{ chNoInput }

workflow peakCallingFlow {

  take:
  bamsChip
  design
  effgsize
  chrsize
  gtf
  fasta
  
  main:
  chVersions = Channel.empty()

  samtoolsFlagstat(
    bamsChip.map{it -> [it[0], it[1]]}
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  // all pairs of samples + no Input
  bamsChip 
    .combine(chNoInput.concat(bamsChip))
    .set { bamsPairChip }  

  //[group, peakType, prefix1, bam1, bai1, prefix2, bam2, bai2]
  design
    .combine(bamsPairChip)
    .filter { it[0] == it[4] && it[1] == it[7] }
    .map { it ->  it[2..-1] }
    .set { chBamCallPeaks }  


  /*********************************
   * Macs2 - sharp mode
   */ 

  chBamCallPeaks
    .filter { it[1] == 'sharp' }
    .map{ it -> it[2..7] }
    .set { chBamMacs2Sharp }

  macs2Sharp(
    chBamMacs2Sharp,
    effgsize.first(),
    chPeakCountHeader.collect()
  )
  chVersions = chVersions.mix(macs2Sharp.out.versions)

  /*******************************
   * Macs2  - Broad
   */

  chBamCallPeaks
    .filter { it[1] == 'broad' }
    .map{ it -> it[2..7] }
    .set { chBamMacs2Broad }

  macs2Broad(
   chBamMacs2Broad,
   effgsize.first(),
   chPeakCountHeader.collect()
  )
  chVersions = chVersions.mix(macs2Broad.out.versions)
     
  /********************************
   * EPIC2 - very broad
   */       

  chBamCallPeaks
    .filter { it[1] == 'very-broad' }
    .map{ it -> it[2..7] }
    .set { chBamEpic }

  epic2(
    chBamEpic,
    effgsize.first(),
    chrsize.collect(),
    chPeakCountHeader.collect()
  )
  chVersions = chVersions.mix(epic2.out.versions)

  // Join the results of all peaks callers
  macs2Sharp.out.peaks
    .mix(macs2Broad.out.peaks, epic2.out.peaks)
    .set{ chPeaks }

  /********************************
   * FRIP
   */

  //[prefix, bam, stats, peaks]
  chBamCallPeaks
    .map{ it -> it[2..7] }
    .combine(samtoolsFlagstat.out.stats, by:0)
    .join(chPeaks)
    .map{it -> [it[0], it[1], it[6], it[7]]}
    .set{ chFrip }

  frip(
    chFrip,
    chFripScoreHeader.collect()
  )
  chVersions = chVersions.mix(frip.out.versions)

  /********************************
   * Peaks Annotation
   */       

  annotPeaksHomer(
    chPeaks,
    gtf.collect(),
    fasta.collect()
  )

  /*******************************
   * QC
   */
       
  peakQC(
    chPeaks.map{it->it[1]}.collect(),
    annotPeaksHomer.out.output.map{it->it[1]}.collect(),
    chPeakAnnotationHeader
  )
  chVersions = chVersions.mix(peakQC.out.versions)

  /*
   * Irreproducible Discovery Rate
   */
       
//       chPeaks 
//         .map { it -> [it[0],it[4]] }
//         .groupTuple()
//         .dump (tag:'rep')
//         .set{ chPeaksPerGroup }
 
//       IDR(
//         chPeaksPerGroup 
//       )

  emit:
    versions = chVersions    
    peaksCountsMqc = macs2Sharp.out.mqc.concat(macs2Broad.out.mqc, epic2.out.mqc)
    peaksOutput = macs2Sharp.out.outputXls.concat(macs2Broad.out.outputXls, epic2.out.output)
    peaksQCMqc = peakQC.out.mqc
    fripResults = frip.out.output
}


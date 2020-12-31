/***********************
 * Peak calling
 */

/* 
 * include requires tasks 
 */
include { sharpMACS2 } from '../processes/sharpMACS2' 
include { broadMACS2 } from '../processes/broadMACS2'
include { veryBroadEpic2 } from '../processes/veryBroadEpic2'
include { peakAnnoHomer } from '../processes/peakAnnoHomer'
include { peakQC } from '../processes/peakQC'
// include { IDR } from '../processes/IDR'

Channel
  .fromPath("$baseDir/assets/peak_count_header.txt")
  .set { chPeakCountHeader }

Channel
  .fromPath("$baseDir/assets/frip_score_header.txt")
  .set { chFripScoreHeader }

Channel
  .fromPath("$baseDir/assets/peak_annotation_header.txt")
  .set{ chPeakAnnotationHeader }

// Chromosome size file
if ( params.chrsize ){
  Channel
    .fromPath(params.chrsize, checkIfExists: true)
    .set{chChromSize}
}
else{
  exit 1, "Chromosome size file not found: ${params.chrsize}"
}

// Annotations
if (params.gtf) {
  Channel
    .fromPath(params.gtf, checkIfExists: true)
    .set{chGtfHomer}
}
else {
  exit 1, "GTF annotation file not specified!"
}

// Fasta file
if ( params.fasta ){
  Channel
    .fromPath(params.fasta, checkIfExists: true)
    .set{chFastaHomer}
}
else{
  exit 1, "Fasta file not found: ${params.fasta}"
}

workflow peakCallingFlow {
    // required inputs
    take:
     chBamsChip 
     chDesignControl
     chNoInput
     chFlagstatMacs 
     
    // workflow implementation
    main:

      /*
       * Prepare channels
       */
      if (params.design){
        chBamsChip 
         .join(chFlagstatMacs)
         .combine(chNoInput.concat(chBamsChip))
         .set { chBamsChip }  

        chDesignControl
         .combine(chBamsChip)
         .filter { it[0] == it[5] && it[1] == it[8] }
         .map { it ->  it[2..-1] }
         .dump(tag:'peakCall')
         .set { chGroupBamMacs }  
       }else{
         chGroupBamMacs=Channel.empty()
       }

       
 
      /*
       * MACS2 - sharp mode
       */
      
      // chGroupBamMacs.filter { it[2] == 'sharp' }.view()
      sharpMACS2(
	chGroupBamMacs.filter { it[2] == 'sharp' },
	chPeakCountHeader.collect(),
	chFripScoreHeader.collect()
      )
      
      /*
       * MACS2  - Broad
       */
      broadMACS2(
        chGroupBamMacs.filter { it[2] == 'broad' },
        chPeakCountHeader.collect(),
        chFripScoreHeader.collect()
      )
      
      /*
       * EPIC2 - very broad
       */

      veryBroadEpic2(
        chGroupBamMacs.filter { it[2] == 'very-broad' },
        chPeakCountHeader.collect(),
        chFripScoreHeader.collect(),
	chChromSize.collect()
      ) 

      // Join the results of all peaks callers
      sharpMACS2.out.peaksMacsSharp  
        .mix(broadMACS2.out.peaksMacsBroad, veryBroadEpic2.out.peaksEpic)
        .set{ chPeaks }

      /*
       * Peaks Annotation
       */

      peakAnnoHomer(
         chPeaks,
         chGtfHomer.collect(),
         chFastaHomer.collect()
       )

      /*
       * Peak calling & annotation QC
       */
       peakQC(
         chPeaks.collect{ it[-1] },
         peakAnnoHomer.out.homerMqc.collect(),
         chPeakAnnotationHeader
       )


    emit:

      version = broadMACS2.out.version // channel: [ path("v_macs2.txt") ]
}


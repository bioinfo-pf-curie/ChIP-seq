/* 
 * STAR Workflow
 * Able to align both on reference and spike genomes
 */

include { starAlign as starAlignRef } from '../process/starAlign'
include { starAlign as starAlignSpike } from '../process/starAlign'
include { compareBams } from '../process/compareBams'
include { samtoolsSort } from '../../common/process/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtoolsFlagstat'

workflow mappingStarFlow {

  take:
  reads 
  indexRef
  genome
  indexSpike
  spike

  main:
  
  chVersions = Channel.empty()
  chGtf = Channel.empty()

  starAlignRef(
    reads,
    indexRef.collect(),
    chGtf.collect().ifEmpty([]),
    Channel.value(false)
  )
  chVersions = chVersions.mix(starAlignRef.out.versions)

  if (spike){

    // Align on spike genome
    starAlignSpike(
      reads,
      indexSpike.collect(),
      chGtf.collect().ifEmpty([]),
      Channel.value(false)
    )

    // Compare reference/spike mapping
    compareBams(starAlignRef.out.bam.join(starAlignSpike.out.bam), genome, spike)

    chAllBams = compareBams.out.refBam
      .concat(compareBams.out.spikeBam)
  }else{
    chAllBams = starAlignRef.out.bam
  }

  // Process all bams
  samtoolsSort(
    chAllBams
  )
  chVersions = chVersions.mix(samtoolsSort.out.versions)

  samtoolsIndex(
    samtoolsSort.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  samtoolsFlagstat(
    samtoolsSort.out.bam
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  emit:
  bam = samtoolsSort.out.bam
  bai = samtoolsIndex.out.bai
  flagstat = samtoolsFlagstat.out.stats
  logs = starAlignRef.out.logs
  versions = chVersions
}

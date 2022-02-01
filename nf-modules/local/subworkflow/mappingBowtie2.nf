/* 
 * Bowtie2 Worflow
 * Able to align both on reference and spike genomes
 */

include { bowtie2 as bowtie2Ref } from '../process/bowtie2'
include { bowtie2 as bowtie2Spike } from '../process/bowtie2'
include { compareBams } from '../process/compareBams'
include { samtoolsSort } from '../../common/process/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtoolsFlagstat'

workflow mappingBowtie2Flow {

  take:
  reads
  indexRef
  genome
  indexSpike
  spike

  main:
  chVersions = Channel.empty()

  bowtie2Ref(
    reads,
    indexRef.collect()
  )
  chVersions = chVersions.mix(bowtie2Ref.out.versions)

  if (spike){

    // Align on spike genome
    bowtie2Spike(
      reads,
      indexSpike.collect()
    )

    // Compare reference/spike mapping
    compareBams(bowtie2Ref.out.bam.join(bowtie2Spike.out.bam), genome, spike)

    chAllBams = compareBams.out.refBam
      .concat(compareBams.out.spikeBam)
  }else{
    chAllBams = bowtie2Ref.out.bam
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
  logs = bowtie2Ref.out.logs
  flagstat = samtoolsFlagstat.out.stats
  versions = chVersions
}

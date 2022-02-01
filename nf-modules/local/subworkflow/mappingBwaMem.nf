/* 
 * BWA Worflow
 * Able to align both on reference and spike genomes
 */

include { bwaMem as bwaMemRef } from '../process/bwaMem'
include { bwaMem as bwaMemSpike } from '../process/bwaMem'
include { compareBams } from '../process/compareBams'
include { samtoolsSort as samtoolsSortRef } from '../../common/process/samtoolsSort'
include { samtoolsSort as samtoolsSortSpike } from '../../common/process/samtoolsSort'
include { samtoolsIndex as samtoolsIndexRef } from '../../common/process/samtoolsIndex'
include { samtoolsIndex as samtoolsIndexSpike } from '../../common/process/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtoolsFlagstat'

workflow mappingBwaMemFlow {

  take:
  reads
  indexRef
  genome
  indexSpike
  spike

  main:
  chVersions = Channel.empty()
  

  // Align on reference genome
  bwaMemRef(
    reads,
    indexRef.collect()
  )
  chVersions = chVersions.mix(bwaMemRef.out.versions)

  if (spike){
    // Align on spike genome
    bwaMemSpike(
      reads,
      indexSpike.collect()
    )

    // Compare reference/spike mapping
    compareBams(bwaMemRef.out.bam.join(bwaMemSpike.out.bam), genome, spike)

    chRefBam = compareBams.out.refBam
    chSpikeBam = compareBams.out.spikeBam
    chCompareBamsMqc = compareBams.out.mqc
    chBwaMemSpikeLogs = bwaMemSpike.out.logs
  }else{
    chRefBam = bwaMemRef.out.bam
    chSpikeBam = Channel.empty()
    chCompareBamsMqc = Channel.empty()
    chBwaMemSpikeLogs = Channel.empty()
  }
 
  // Process reference bams
  samtoolsSortRef(
    chRefBam
  )
  chVersions = chVersions.mix(samtoolsSortRef.out.versions)

  samtoolsIndexRef(
    samtoolsSortRef.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndexRef.out.versions)

  samtoolsFlagstat(
   samtoolsSortRef.out.bam
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  // Process spike bams
  if (spike){
    samtoolsSortSpike(
      chSpikeBam
    )

    samtoolsIndexSpike(
      samtoolsSortSpike.out.bam
    )
    chSpikeBamOutput=samtoolsSortSpike.out.bam.join(samtoolsIndexSpike.out.bai)
  }else{
    chSpikeBamOutput=Channel.empty()
  }

  emit:
  bam = samtoolsSortRef.out.bam.join(samtoolsIndexRef.out.bai)
  logs = bwaMemRef.out.logs
  flagstat = samtoolsFlagstat.out.stats
  spikeBam = chSpikeBamOutput
  spikeLogs = chBwaMemSpikeLogs
  compareBamsMqc = chCompareBamsMqc
  versions = chVersions
}

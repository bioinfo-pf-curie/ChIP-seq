/* 
 * Mapping Worflow
 * Able to align both on reference and spike genomes
 */


if (params.aligner == "bwa-mem"){
  include { bwaMem as mapping } from '../../common/process/bwa/bwaMem'
  include { bwaMem as mappingSpike } from '../../common/process/bwa/bwaMem'
}else if (params.aligner == "bowtie2"){
  include { bowtie2 as mapping } from '../../common/process/bowtie2/bowtie2'
  include { bowtie2 as mappingSpike } from '../../common/process/bowtie2/bowtie2'
}else if (params.aligner == "star"){
  include { starAlign as mapping } from '../../common/process/star/starAlign'
  include { starAlign as mappingSpike } from '../../common/process/star/starAlign'
}

include { compareBams } from '../../local/process/compareBams'
include { samtoolsSort as samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsSort as samtoolsSortSpike } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex as samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsIndex as samtoolsIndexSpike } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'

workflow mappingFlow {

  take:
  reads
  indexRef
  indexSpike

  main:
  chVersions = Channel.empty()
  
  // Align on reference genome
  if (params.aligner == "star"){
    mapping(
      reads,
      indexRef.collect(),
      Channel.empty().collect().ifEmpty([])
    )
  }else{
    mapping(
      reads,
      indexRef.collect()
    )
  }
  chVersions = chVersions.mix(mapping.out.versions)

  if (params.spike){
    // Align on spike genome
    if (params.aligner == "star"){
      mappingSpike(
        reads,
        indexSpike.collect(),
        Channel.empty().collect().ifEmpty([])
      )
    }else{
      mappingSpike(
        reads,
        indexSpike.collect()
      )
    }

    // Compare reference/spike mapping
    compareBams(mapping.out.bam.join(mappingSpike.out.bam), params.genome, params.spike)
    chBam = compareBams.out.refBam
    chSpikeBam = compareBams.out.spikeBam
    chCompareBamsMqc = compareBams.out.mqc
    chSpikeLogs = mappingSpike.out.logs
  }else{
    chBam = mapping.out.bam
    chSpikeBam = Channel.empty()
    chCompareBamsMqc = Channel.empty()
    chSpikeLogs = Channel.empty()
  }
 
  // Add genome information
  chBam = chBam.map{ meta, bam ->
    def newMeta = [ id: meta.id, name: meta.name, singleEnd: meta.singleEnd, genome: params.genome ]
    [newMeta, bam]
  }

  chSpikeBam = chSpikeBam.map{ meta, bam ->
    def newMeta = [ id: meta.id, name: meta.name, singleEnd: meta.singleEnd, genome: params.spike ]
    [newMeta, bam]
  }

  // Process reference bams
  samtoolsSort(
    chBam
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

  // Process spike bams
  if (params.spike){
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
  bam = samtoolsSort.out.bam.join(samtoolsIndex.out.bai)
  logs = mapping.out.logs
  flagstat = samtoolsFlagstat.out.stats
  spikeBam = chSpikeBamOutput
  spikeLogs = chSpikeLogs
  compareBamsMqc = chCompareBamsMqc
  versions = chVersions
}

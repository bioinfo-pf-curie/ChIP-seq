/* 
 * Load BAM files
 */

include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'

workflow loadBamFlow {

  take:
  bam

  main:
  chVersions = Channel.empty()

  // Add genome information
  bam = bam.map{ meta, bam ->
    def newMeta = [ id: meta.id, name: meta.name, singleEnd: meta.singleEnd, genome: params.genome ]
    [newMeta, bam]
  }
  
  // Process reference bams
  samtoolsSort(
    bam
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
  bam = samtoolsSort.out.bam.join(samtoolsIndex.out.bai)
  flagstat = samtoolsFlagstat.out.stats
  versions = chVersions
}

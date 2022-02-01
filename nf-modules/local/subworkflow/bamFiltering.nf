/* 
 * Filter BAMs file
 */

include { markDuplicates } from '../../common/process/markDuplicates' 
include { samtoolsFilter } from '../../common/process/samtoolsFilter'
include { samtoolsFlagstat } from '../../common/process/samtoolsFlagstat'
include { samtoolsIndex } from '../../common/process/samtoolsIndex'

workflow bamFilteringFlow {

    take:
    bams // [prefix, bam, bai]

    main:
    chVersions = Channel.empty()

    // Remove duplicates
    markDuplicates(
      bams
    )
    chVersions = chVersions.mix(markDuplicates.out.versions)

    // Filter reads
    samtoolsFilter(
      markDuplicates.out.bam
    )
    chVersions = chVersions.mix(samtoolsFilter.out.versions)

    // index
    samtoolsIndex(
      samtoolsFilter.out.bam
    )

    // flagstat
    samtoolsFlagstat(
      samtoolsFilter.out.bam
    )

    // idxstats
    // stats

    emit:
    bam  = samtoolsFilter.out.bam.join(samtoolsIndex.out.bai)
    flagstat  = samtoolsFlagstat.out.stats
    markdupMetrics = markDuplicates.out.metrics
    versions = chVersions
}


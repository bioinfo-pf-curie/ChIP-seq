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
    if (!params.skipFiltering){
      samtoolsFilter(
        markDuplicates.out.bam
      )
      chVersions = chVersions.mix(samtoolsFilter.out.versions)
      chBam = samtoolsFilter.out.bam
    }else{
      chBam = markDuplicates.out.bam
    }

    // index
    samtoolsIndex(
      chBam
    )
 
    // flagstat
    samtoolsFlagstat(
      chBam
    )
    

    emit:
    bam  = chBam.join(samtoolsIndex.out.bai)
    flagstat  = samtoolsFlagstat.out.stats
    markdupMetrics = markDuplicates.out.metrics
    versions = chVersions
}


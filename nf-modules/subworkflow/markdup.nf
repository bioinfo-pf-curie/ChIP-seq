/* 
 * Alignment on reference genome 
 */

/* 
 * include requires tasks 
 */
include { markDuplicates } from '../processes/markDuplicates' 
include { preseq } from '../processes/preseq'
include { bamFiltering } from '../processes/bamFiltering'

workflow markdupFlow {
    // required inputs
    take:
      chSortBams
    // workflow implementation
    main:

      markDuplicates(chSortBams)
      preseq(markDuplicates.out.bams.filter{ it[0][-5..-1] != "spike" }) 
      bamFiltering(markDuplicates.out.bams)

    // emit:
    chMarkedPicstats = markDuplicates.out.picstats   // channel: [ path "*metrics.txt" ]
    chPicardVersion = markDuplicates.out.version     // channel: [ path("v_picard.txt") ]
    chPreseqStats = preseq.out.stats                 // channel: [ path("*.ccurve.txt") ]
    chPreseqVersion = preseq.out.version             // channel: [ path("v_preseq.txt") ]
    chFilteredBams  = bamFiltering.out.filteredBams  // channel: [ val(prefix), path("*filtered.{bam,bam.bai}") ]
    chFilteredFlagstat  = bamFiltering.out.filteredFlagstat  // channel: [ val(prefix), path("*filtered.flagstat") ]
    chSamtoolsVersionBamFiltering  = bamFiltering.out.version  // channel: [ path("v_samtools.txt") ]
}


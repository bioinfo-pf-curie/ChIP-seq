/* 
 * Alignment on reference genome 
 */

/* 
 * include requires tasks 
 */
include { markDuplicates } from '../process/markDuplicates' 
include { preseq } from '../process/preseq'
include { bamFiltering } from '../process/bamFiltering'

workflow markdupFlow {
    // required inputs
    take:
      chSortBams
    // workflow implementation
    main:

      markDuplicates(chSortBams)
      bamFiltering(markDuplicates.out.bams)
      preseq(markDuplicates.out.bams.filter{ it[0][-5..-1] != "spike" }) 

    emit:
      chFilteredBams  = bamFiltering.out.filteredBams  // channel: [ val(prefix), path("*filtered.{bam,bam.bai}") ]
      chFilteredFlagstat  = bamFiltering.out.filteredFlagstat  // channel: [ val(prefix), path("*filtered.flagstat") ]
      chMarkedPicstats = markDuplicates.out.picstats   // channel: [ path "*metrics.txt" ]
      chPicardVersion = markDuplicates.out.version     // channel: [ path("v_picard.txt") ]
      chPreseqStats = preseq.out.stats                 // channel: [ path("*.ccurve.txt") ]
      chPreseqVersion = preseq.out.version             // channel: [ path("v_preseq.txt") ]
      chSamtoolsVersionBamFiltering  = bamFiltering.out.version  // channel: [ path("v_samtools.txt") ]
}


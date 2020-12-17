/* 
 * Alignment on reference genome 
 */

/* 
 * include requires tasks 
 */
include { markDuplicates } from '../processes/markDuplicates' 
include { preseq } from '../processes/preseq'

workflow markdupFlow {
    // required inputs
    take:
      chSortBams
    // workflow implementation
    main:

      markDuplicates(chSortBams)
      preseq(markDuplicates.out.bams.filter{ it[0][-5..-1] != "spike" }) 

    // emit:
    //  bam // channel: [ val(sample), path("*.bam") ]
    //  mqc // channel: [ path("*.log") ]
}


/* 
 * Alignment on reference genome 
 */

/* 
 * include requires tasks 
 */
include { bwaMem } from '../processes/bwaMem' 
include { bowtie2 } from '../processes/bowtie2'
include { star } from '../processes/star'

workflow mappingFlow {
    // required inputs
    take:
      rawReads 
      chBwaIndex
      chBt2Index
      chStarIndex
    // workflow implementation
    main:

      if (params.aligner == "bowtie2"){
        bowtie2(rawReads.combine(chBt2Index))
        bam = bowtie2.out.bam
        mqc = bowtie2.out.mqc
      } else if (params.aligner == "bwa-mem"){
        bwaMem(rawReads.combine(chBwaIndex))
        bam = bwaMem.out.bam
        mqc = bwaMem.out.mqc
      } else if (params.aligner == "star"){
        star(rawReads.combine(chStarIndex))
        bam = star.out.bam
        mqc = star.out.mqc
      }

    emit:
      bam // channel: [ val(sample), path("*.bam") ]
      mqc // channel: [ path("*.log") ]
}


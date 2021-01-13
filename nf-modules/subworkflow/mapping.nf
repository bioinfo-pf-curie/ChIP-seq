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
      chBwaVersion = Channel.empty()
      chBowtie2Version = Channel.empty()
      chStarVersion = Channel.empty()
      if (params.aligner == "bowtie2"){
        bowtie2(rawReads.combine(chBt2Index))
        bam = bowtie2.out.bam
        mqc = bowtie2.out.mqc
        chBowtie2Version = bowtie2.out.version
      } else if (params.aligner == "bwa-mem"){
        bwaMem(rawReads.combine(chBwaIndex))
        bam = bwaMem.out.bam
        mqc = bwaMem.out.mqc
        chBwaVersion = bwaMem.out.version
      } else if (params.aligner == "star"){
        star(rawReads.combine(chStarIndex))
        bam = star.out.bam
        mqc = star.out.mqc
        chStarVersion = star.out.version
      }

    emit:
      bam // channel: [ val(sample), path("*.bam") ]
      mqc // channel: [ path("*.log") ]
      chBwaVersion // channel: [ path("v_bwa.txt") ]
      chBowtie2Version // channel: [ path("v_bowtie2.txt") ]
      chStarVersion  // channel: [ path("v_star.txt") ]
}


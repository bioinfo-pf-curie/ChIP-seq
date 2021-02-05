/* 
 * Alignment on reference genome 
 */

/* 
 * include requires tasks 
 */
include { bwaMem } from '../process/bwaMem' 
include { bowtie2 } from '../process/bowtie2'
include { star } from '../process/star'

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
      chAlignReads = Channel.empty()
      chMappingMqc = Channel.empty()
      if (params.aligner == "bowtie2"){
        bowtie2(rawReads.combine(chBt2Index))
        chAlignReads = bowtie2.out.bam
        chMappingMqc = bowtie2.out.mqc
        chBowtie2Version = bowtie2.out.version
      } else if (params.aligner == "bwa-mem"){
        bwaMem(rawReads.combine(chBwaIndex))
        chAlignReads = bwaMem.out.bam
        chMappingMqc = bwaMem.out.mqc
        chBwaVersion = bwaMem.out.version
      } else if (params.aligner == "star"){
        star(rawReads.combine(chStarIndex))
        chAlignReads = star.out.bam
        chMappingMqc = star.out.mqc
        chStarVersion = star.out.version
      }

    emit:
      chAlignReads     // channel: [ val(sample), path("*.bam") ]
      chMappingMqc     // channel: [ path("*.log") ]
      chBwaVersion     // channel: [ path("v_bwa.txt") ]
      chBowtie2Version // channel: [ path("v_bowtie2.txt") ]
      chStarVersion    // channel: [ path("v_star.txt") ]
}


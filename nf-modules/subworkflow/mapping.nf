// genomeRef = params.genome
/* 
 * include requires tasks 
 */
include { bwaMem } from '../processes/bwaMem' 
//include { bwaMem } from '../processes/bwaMem' addParams( genomeRef: genomeRef)
include { bowtie2 } from '../processes/bowtie2'
include { star } from '../processes/star'

/* 
 * Alignment on reference genome 
 */
workflow mappingFlow {
    // required inputs
    take:
      rawReads 
      chBwaIndex
      chBt2Index
      chStarIndex
    // workflow implementation
    main:
      chAlignReads = Channel.empty()
      chMappingMqc = Channel.empty()
      bwaMem(rawReads.combine(chBwaIndex))
      bowtie2(rawReads.combine(chBt2Index))
      star(rawReads.combine(chStarIndex))

      if (params.aligner == "bowtie2"){
        bam = bowtie2.out.bam
        mqc = bowtie2.out.mqc
      } else if (params.aligner == "bwa-mem"){
        bam = bwaMem.out.bam
        mqc = bwaMem.out.mqc
      } else if (params.aligner == "star"){
        bam = star.out.bam
        mqc = star.out.mqc
      }

    emit:
      bam // channel: [ val(sample), path("*.bam") ]
      mqc // channel: [ path("*.log") ]
}


/*
 * Alignment on reference genome
 */

/* BWA-MEM */
process bwaMem{
  tag "${sample} on ${genomeName}"
  label 'bwa'
  label 'highCpu'
  label 'highMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
             saveAs: {filename ->
             if (filename.indexOf(".log") > 0) "logs/$filename"
             else if (params.saveAlignedIntermediates) filename}

  when:
  params.aligner == "bwa-mem" && !params.inputBam

  input:
  val genomeRef
  tuple val(sample), path(reads), path(index), val(genomeBase), val(genomeName)

  output:
  tuple val(sample), path("*.bam"), emit: bam 
  path "*.log"                    , emit: mqc 
  path "v_bwa.txt"                , emit: version

  script:
  prefix = genomeName == genomeRef ? sample : sample + '_spike'
  opts = params.bwaOpts
  """
  echo \$(bwa 2>&1) &> v_bwa.txt
  bwa mem -t ${task.cpus} \
           ${index}/${genomeBase} \
          ${opts} \
          $reads | samtools view -bS - > ${prefix}.bam
  getBWAstats.sh -i ${prefix}.bam -p ${task.cpus} > ${prefix}_bwa.log
  """
}


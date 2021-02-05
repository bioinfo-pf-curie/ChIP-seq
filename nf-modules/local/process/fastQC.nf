/*
 * FastQC
 */
process fastQC{
  tag "${prefix}"
  label 'fastqc'
  label 'lowCpu'
  label 'minMem'
  publishDir "${params.outDir}/fastqc", mode: 'copy'

  when:
  !params.skipFastqc && !params.inputBam

  input:
  tuple val(prefix), path(reads) 

  output:
  path "*_fastqc.{zip,html}", emit: mqc 
  path "v_fastqc.txt"       , emit: version

  script:
  """
  fastqc --version &> v_fastqc.txt
  fastqc -q $reads -t ${task.cpus}
  """
}


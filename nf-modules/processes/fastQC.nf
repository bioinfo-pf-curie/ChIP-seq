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
  set val(prefix), file(reads) 

  output:
  file "*_fastqc.{zip,html}" 
  file("v_fastqc.txt")

  script:
  """
  fastqc --version &> v_fastqc.txt
  fastqc -q $reads -t ${task.cpus}
  """
}


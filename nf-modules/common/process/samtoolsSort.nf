/*
 * Samtools - Sort
 */

process samtoolsSort {
  tag "$prefix"
  label 'samtools'
  label 'medCpu'
  label 'medMem'
 
  input:
  tuple val(prefix), path (bam)

  output:
  tuple val(prefix), path ("*_sorted.bam"), emit: bam
  path("versions.txt") , emit: versions

  script:
  def args = task.ext.args ?: ''
  """
  echo \$(samtools --version | head -1 ) > versions.txt
  samtools sort \\
    ${args} \\
    -@  ${task.cpus}  \\
    -o ${bam.baseName}_sorted.bam  \\
    ${bam}
  """
}
    
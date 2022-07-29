/*
 * Samtools - Sort
 */

process samtoolsSort {
  tag "${meta.id}"
  label 'samtools'
  label 'medCpu'
  label 'medMem'
 
  input:
  tuple val(meta), path (bam)

  output:
  tuple val(meta), path ("*_sorted.bam"), emit: bam
  path("versions.txt") , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${bam.baseName}"
  """
  echo \$(samtools --version | head -1 ) > versions.txt
  samtools sort \\
    ${args} \\
    -@  ${task.cpus}  \\
    -o ${prefix}_sorted.bam  \\
    ${bam}
  """
}
    
/*
 * BAM Filtering based on samtools
 */

process samtoolsFilter {
  tag "${meta.id}"
  label 'samtools'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("*filtered.bam"), emit: bam
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${bam.baseName}"
  """
  echo \$(samtools --version | head -1) > versions.txt
  samtools view \\
    ${args} \\
    -b ${bam} > ${prefix}_filtered.bam
  """
}

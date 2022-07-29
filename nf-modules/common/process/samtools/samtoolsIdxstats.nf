/*
 * Samtools - Idxstats
 */

process samtoolsIdxstats {
  tag "${meta.id}"
  label 'samtools'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(meta), path (bam)

  output:
  tuple val(meta), path("*idxstats"), emit: stats
  path("versions.txt") , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(samtools --version | head -1) > versions.txt
  samtools idxstats ${bam} > ${prefix}.idxstats
  """
}

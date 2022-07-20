/*
 * Get fragment sizes with Picard
 */

process collectInsertSizeMetrics {
  tag "${meta.id}"
  label 'picard'
  label 'lowCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  path("*insert_size*.{pdf,txt}"), emit: results
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(picard CollectInsertSizeMetrics --version 2>&1 | sed -e 's/Version:/picard /') > versions.txt
  picard CollectInsertSizeMetrics \\
      I=${bam} \\
      O=${prefix}_insert_size_metrics.txt \\
      H=${prefix}_insert_size_histogram.pdf \\
      ${args}
  """
}


/*
 * Get fragment sizes
 */

process getFragmentSize {
  tag "${prefix}"
  label 'picard'
  label 'lowCpu'
  label 'medMem'

  input:
  tuple val(prefix), path(bam), path(bai)

  output:
  path("*.{pdf,txt}"), emit: fragmentsSize
  path("versions.txt"), emit: versions

  script:
  """
  echo \$(picard CollectInsertSizeMetrics --version 2>&1 | sed -e 's/Version:/picard /') > versions.txt
  picard CollectInsertSizeMetrics \
      I=${bam} \
      O=${prefix}_insert_size_metrics.txt \
      H=${prefix}_insert_size_histogram.pdf \
      VALIDATION_STRINGENCY=LENIENT \
      M=0.5
  """
}


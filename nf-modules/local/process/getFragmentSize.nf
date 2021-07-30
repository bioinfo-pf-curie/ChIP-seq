/*
 * Get fragment sizes
 */

process getFragmentSize {
  tag "${prefix}"
  label 'picard'
  label 'lowCpu'
  label 'medMem'

  publishDir path: "${params.outDir}/fragSize/", mode: "copy"

  when:
  !params.singleEnd

  input:
  tuple val(prefix), path(filteredBam) 

  output:
  path("*.{pdf,txt}"),   emit: fragmentsSize

  script:
  """
  picard CollectInsertSizeMetrics \
      I=${filteredBam[0]} \
      O=${prefix}_insert_size_metrics.txt \
      H=${prefix}_insert_size_histogram.pdf \
      VALIDATION_STRINGENCY=LENIENT \
      M=0.5
  """
}


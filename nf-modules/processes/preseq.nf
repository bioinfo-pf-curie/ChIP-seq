/*
 * Preseq (before alignment filtering and only on ref mapped reads)
 */

process preseq {
  tag "${prefix}"
  label 'preseq'
  label 'lowCpu'
  label 'minMem'
  publishDir "${params.outDir}/preseq", mode: 'copy'

  when:
  !params.skipPreseq

  input:
  tuple val(prefix), path(bam)

  output:
  path "*.ccurve.txt" 
  path("v_preseq.txt")

  script:
  defectMode = params.preseqDefect ? '-D' : ''
  """
  preseq &> v_preseq.txt
  preseq lc_extrap -v $defectMode -output ${prefix}.ccurve.txt -bam ${bam[0]} -e 200e+06
  """
}


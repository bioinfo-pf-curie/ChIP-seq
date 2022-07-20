/*
 * Feature counts on BED file
 */

process bed2saf{
  label 'unix'
  label 'minCpu'
  label 'minMem'

  input:
  path(bed)

  output:
  path("*.saf"), emit: saf

  script:
  def prefix = task.ext.prefix ?: "${bed.baseName}"
  """
  echo "GeneID\tChr\tStart\tEnd\tStrand" > ${prefix}.saf
  awk '{OFS="\t";print \$4,\$1,\$2,\$3,\$6}' ${bed} >> ${prefix}.saf
  """
}

/**************************************
 * Feature counts on BED file
 */

process featureCounts{
  label 'featureCounts'
  label 'medCpu'
  label 'lowMem'
  tag("${annot}") 

  input:
  path(bams)
  each path(annot)

  output:
  path("*csv"), emit: counts
  path("*summary"), emit: logs
  path("versions.txt"), emit: versions 

  script:
  prefix = annot.toString() - ~/(\.bed)?$/
  def args = task.ext.args ?: ''
  """
  echo \$(featureCounts -v 2>&1 | sed '/^\$/d') > versions.txt
  awk '{OFS="\t";print \$4,\$1,\$2,\$3,\$6}' ${annot} > ${prefix}.saf
  featureCounts -a ${prefix}.saf -F SAF \\
                -o allchip_counts_${prefix}.csv \\
                -T ${task.cpus} \\
                -s 0 \\
		${args} \\
                -O ${bams} 2> featureCounts_${prefix}.log
  """

}

/**************************************
 * Feature counts
 */

process featureCounts{
  tag "${bed}"
  label 'featureCounts'
  label 'medCpu'
  label 'lowMem'
  publishDir "${params.outDir}/featCounts/", mode: "copy"

  when:
  !params.skipFeatCounts

  input:
  path(bams)
  each path(annot)

  output:
  path("*csv")               , emit: featCounts
  path("*summary")           , emit: featCountsMqc
  path("v_featurecounts.txt"), emit: version 

  script:
  prefix = annot.toString() - ~/(\.bed)?$/
  paramsPairedEnd = params.singleEnd ? '' : '-p -C -B'
  """
  featureCounts -v &> v_featurecounts.txt
  awk '{OFS="\t";print \$4,\$1,\$2,\$3,\$6}' ${annot} > ${prefix}.saf
  featureCounts -a ${prefix}.saf -F SAF \\
                -o allchip_counts_${prefix}.csv \\
                -T ${task.cpus} \\
                -s 0 ${paramsPairedEnd} \\
                -O ${bams} 2> featureCounts_${prefix}.log
  """

}

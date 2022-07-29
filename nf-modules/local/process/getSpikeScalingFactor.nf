/*
 * Calculate spike-in scaling factors
 */

process getSpikeScalingFactor {
  label 'r'
  label 'minCpu'
  label 'medMem'

  input:
  path(tab)

  output:
  path("*.sf"), emit: tabSF
  path("versions.txt"), emit: versions

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') >> versions.txt
  getDESeqSF.r ${tab}
  """
}

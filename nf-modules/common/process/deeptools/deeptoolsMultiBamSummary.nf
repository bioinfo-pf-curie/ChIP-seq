/*
 * Deeptools - multiBamSummary
 */

process deeptoolsMultiBamSummary {
  label 'deeptools'
  label 'medCpu'
  label 'medMem'
  
  input:
  path(bams)
  path(bai)

  output:
  path("*readCounts.tab"), emit: output
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "multiBam"
  """
  echo \$(deeptools --version ) > versions.txt
  multiBamSummary bins \\
                 -b $bams \\
		 ${args} \\
                 -o results.npz \\
                 --outRawCounts ${prefix}_readCounts.tab
  """
}


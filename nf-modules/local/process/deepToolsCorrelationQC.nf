/*
 * Deeptools - correlatation
 */

process deepToolsCorrelationQC{
  label 'deeptools'
  label 'highCpu'
  label 'lowMem'

  when:
  allPrefix.size() >= 2

  input:
  val(allPrefix)
  path(allBams) 
  path(allBai) 
  path(BLbed)

  output:
  path("bams_correlation.{pdf,tab}"), emit: output
  path("*tab"), emit: mqc
  path("versions.txt"), emit: versions

  script:
  allPrefix = allPrefix.toString().replace("[","")
  allPrefix = allPrefix.replace(","," ")
  allPrefix = allPrefix.replace("]","")
  def args = task.ext.args ?: ''
  blacklist = params.blacklist ? "--blackListFileName ${BLbed}" : ""
  """
  echo \$(deeptools --version ) > versions.txt
  multiBamSummary bins -b $allBams \\
                       ${args} \\
		       ${blacklist} \\
                       -o bams_summary.npz \\
                       -p ${task.cpus}

  plotCorrelation -in bams_summary.npz \\
                  -o bams_correlation.pdf \\
                  -c spearman -p heatmap -l $allPrefix \\
                  --outFileCorMatrix bams_correlation.tab
  """
}


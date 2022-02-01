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
  blacklistParams = params.blacklist ? "--blackListFileName ${BLbed}" : ""
  allPrefix = allPrefix.toString().replace("[","")
  allPrefix = allPrefix.replace(","," ")
  allPrefix = allPrefix.replace("]","")
  """
  echo \$(deeptools --version ) > versions.txt
  multiBamSummary bins -b $allBams \\
                       --binSize=50000 \\
                        -o bams_summary.npz \\
                        -p ${task.cpus} \\
                        ${blacklistParams}

  plotCorrelation -in bams_summary.npz \\
                  -o bams_correlation.pdf \\
                  -c spearman -p heatmap -l $allPrefix \\
                  --outFileCorMatrix bams_correlation.tab
  """
}


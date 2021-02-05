process deepToolsCorrelationQC{
  label 'deeptools'
  label 'highCpu'
  label 'lowMem'
  publishDir "${params.outDir}/deepTools/correlationQC", mode: "copy"

  when:
  allPrefix.size() >= 2 && !params.skipDeepTools

  input:
  path(allBams) 
  path(allBai) 
  val (allPrefix)
  path(BLbed)

  output:
  path "bams_correlation.pdf", emit: correl
  path "bams_correlation.tab", emit: correlMqc

  script:
  blacklistParams = params.blacklist ? "--blackListFileName ${BLbed}" : ""
  allPrefix = allPrefix.toString().replace("[","")
  allPrefix = allPrefix.replace(","," ")
  allPrefix = allPrefix.replace("]","")
  """
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


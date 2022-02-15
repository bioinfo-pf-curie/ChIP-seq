/*
 * Deeptools - FingerPrints
 */

process deepToolsFingerprint{
  label 'deeptools'
  label 'highCpu'
  label 'lowMem'

  input:
  val(allPrefix)
  path(allBams) 
  path(allBai) 

  output:
  path("*.{pdf,txt}"), emit: output
  path("*.txt"), emit: mqc
  path("versions.txt"), emit: versions

  script:
  allPrefix = allPrefix.toString().replace("[","")
  allPrefix = allPrefix.replace(","," ")
  allPrefix = allPrefix.replace("]","")
  def args = task.ext.args ?: ''
  """
  echo \$(deeptools --version ) > versions.txt
  plotFingerprint -b $allBams \\
                  -plot bams_fingerprint.pdf \\
                  -p ${task.cpus} \\
                  -l $allPrefix \\
                  ${args} \\
                  --outRawCounts plotFingerprint.raw.txt \\
                  --outQualityMetrics plotFingerprint.qmetrics.txt
  """
}



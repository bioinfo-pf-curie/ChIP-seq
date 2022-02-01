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
  if (params.singleEnd){
    extend = params.fragmentSize > 0 ? "--extendReads ${params.fragmentSize}" : ""
  }else{
    extend = params.noReadExtension ? "" : "--extendReads"
  }
  allPrefix = allPrefix.toString().replace("[","")
  allPrefix = allPrefix.replace(","," ")
  allPrefix = allPrefix.replace("]","")
  """
  echo \$(deeptools --version ) > versions.txt
  plotFingerprint -b $allBams \\
                  -plot bams_fingerprint.pdf \\
                  -p ${task.cpus} \\
                  -l $allPrefix \\
                  ${extend} \\
                  --skipZeros \\
                  --outRawCounts plotFingerprint.raw.txt \\
                  --outQualityMetrics plotFingerprint.qmetrics.txt
  """
}

process deepToolsFingerprint{
  label 'deeptools'
  label 'highCpu'
  label 'lowMem'
  publishDir "${params.outDir}/deepTools/fingerprintQC", mode: "copy"

  when:
  !params.skipDeepTools

  input:
  path(allBams) 
  path(allBai) 
  val (allPrefix)

  output:
  path "bams_fingerprint.pdf", emit: fingerprint
  path "plotFingerprint*"    , emit: fingerprintMqc

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

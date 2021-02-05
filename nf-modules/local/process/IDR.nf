process IDR{
  tag "${group}"
  label 'idr'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/IDR", mode: 'copy'
  
  when:
  allPeaks.toList().size > 1 && !params.skipIDR
  
  input:
  tuple val(group), path(allPeaks)
  
  output:
  path "*idrValues.txt", emit: idr 
  path "*log.txt"      , emit: mqcIdr 
  path("v_idr.txt")    , emit: version
  
  script:
  peaktype = allPeaks[0].toString()
  peaktype = peaktype.substring(peaktype.lastIndexOf(".") + 1)
  """ 
  idr --version &> v_idr.txt
  idr --samples ${allPeaks} \\
      --input-file-type  ${peaktype} \\
      -o ${group}_idrValues.txt \\
      -l ${group}_log.txt \\
      --plot
  """
}

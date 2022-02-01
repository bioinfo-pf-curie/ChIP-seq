/*
 * Macs2 - peak calling
 */

process macs2{
  tag "${chipPrefix} - ${controlPrefix}"
  label 'macs2'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(chipPrefix), path(chipBam), path(chipBai), val(controlPrefix), path(controlBam), path(controlBai)
  val(effGenomeSize)
  path peakCountHeader

  output:
  path("*.xls"), emit: outputXls
  tuple val(chipPrefix), path("*.{narrowPeak,broadPeak}"), emit: peaks
  path("*_mqc.tsv"), emit: mqc
  path("versions.txt"), emit: versions

  script:
  format = params.singleEnd ? "BAM" : "BAMPE"
  ctrl = controlPrefix != 'NO_INPUT' ? "-c ${controlBam}" : ''
  def args = task.ext.args ?: ''
  def outputSuffix = (args.contains('--broad')) ? "broadPeak" : "narrowPeak"
  """
  echo \$(macs2 --version 2>&1) &> versions.txt
  macs2 callpeak \\
    -t ${chipBam} \\
    ${ctrl} \\
    -f $format \\
    -n ${chipPrefix}_${controlPrefix}_macs2 \\
    -g $effGenomeSize \\
    ${args}

  cat ${chipPrefix}_${controlPrefix}_macs2_peaks.${outputSuffix} | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${chipPrefix}", \$1 }' | cat $peakCountHeader - > ${chipPrefix}_${controlPrefix}_macs2_peaks.count_mqc.tsv
  """
}



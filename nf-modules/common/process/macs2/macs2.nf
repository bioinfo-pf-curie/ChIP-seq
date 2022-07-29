/*
 * Macs2 - peak calling
 */

process macs2{
  tag "${prefix}"
  label 'macs2'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai), path(controlBam), path(controlBai)
  val(effGenomeSize)
  path(peakCountHeader)

  output:
  path("*.xls"), emit: outputXls
  tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peaks
  path("*_mqc.tsv"), emit: mqc
  path("versions.txt"), emit: versions

  script:
  format = meta.singleEnd ? "BAM" : "BAMPE"
  ctrl = controlBam ? "-c ${controlBam}" : ''
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  def outputSuffix = (args.contains('--broad')) ? "broadPeak" : "narrowPeak"
  """
  echo \$(macs2 --version 2>&1) &> versions.txt
  macs2 callpeak \\
    ${args} \\
    -t ${bam} \\
    ${ctrl} \\
    -f $format \\
    -n ${prefix}_macs2 \\
    -g $effGenomeSize \\

  cat ${prefix}_macs2_peaks.${outputSuffix} | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${meta.id}", \$1 }' | cat $peakCountHeader - > ${prefix}_macs2_peaks.count_mqc.tsv
  """
}



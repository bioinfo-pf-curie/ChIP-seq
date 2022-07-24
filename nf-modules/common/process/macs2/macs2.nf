/*
 * Macs2 - peak calling
 */

process macs2{
  tag "${meta.id} - ${meta.control}"
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
  ctrl = meta.control != 'NO_INPUT' ? "-c ${controlBam}" : ''
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}_${meta.control}"
  def outputSuffix = (args.contains('--broad')) ? "broadPeak" : "narrowPeak"
  """
  echo \$(macs2 --version 2>&1) &> versions.txt
  macs2 callpeak \\
    -t ${bam} \\
    ${ctrl} \\
    -f $format \\
    -n ${prefix}_macs2 \\
    -g $effGenomeSize \\
    ${args}

  cat ${prefix}_macs2_peaks.${outputSuffix} | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${meta.id}", \$1 }' | cat $peakCountHeader - > ${prefix}_macs2_peaks.count_mqc.tsv
  """
}



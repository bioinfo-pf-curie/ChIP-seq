/*
 * EPIC2 - peak calling for large histone marks
 */

process epic2{
  tag "${chipPrefix} - ${controlPrefix}"
  label 'epic2'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(chipPrefix), path(chipBam), path(chipBai), val(controlPrefix), path(controlBam), path(controlBai)
  val(effGenomeSize)
  path(chromsize)
  path peakCountHeader

  output:
  tuple val(chipPrefix), path("*.broadPeak"), emit: peaks
  path("*.out"), emit: output
  path("*_mqc.tsv"), emit: mqc
  path("versions.txt"), emit: versions

  script:
  ctrl = controlPrefix != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  def args = task.ext.args ?: ''
  """
  echo "epic2 "\$(epic2 --version) &> versions.txt
  epic2 -t ${chipBam} \\
    ${ctrl} \\
    --chromsizes ${chromsize} \\
    --effective-genome-fraction ${effGenomeSize} \\
    -o ${chipPrefix}_${controlPrefix}_epic.out \\
    ${args}

  echo "track type=broadPeak name=\"${chipPrefix}\" description=\"${chipPrefix}_${controlPrefix}\" nextItemButton=on" > ${chipPrefix}_${controlPrefix}_epic2_peaks.broadPeak
  awk -v id="${chipPrefix}" 'NF<10{fc=0;pval=-1;qval=-1}NF==10{fc=\$10;pval=\$4;qval=\$9}NR>1{OFS="\t"; print \$1,\$2,\$3,id"_peak_"NR,\$5,".",fc,pval,qval}' ${chipPrefix}_${controlPrefix}_epic.out >> ${chipPrefix}_${controlPrefix}_epic2_peaks.broadPeak
  cat ${chipPrefix}_${controlPrefix}_epic2_peaks.broadPeak | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${chipPrefix}", \$1 }' | cat $peakCountHeader - > ${chipPrefix}_${controlPrefix}_epic2_peaks.count_mqc.tsv
  """
}




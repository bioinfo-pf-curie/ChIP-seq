/*
 * EPIC2 - peak calling for large histone marks
 */

process epic2{
  tag "${meta.id} - ${meta.control}"
  label 'epic2'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai), path(controlBam), path(controlBai)
  val(effGenomeSize)
  path(chromsize)
  path(peakCountHeader)

  output:
  tuple val(meta), path("*.broadPeak"), emit: peaks
  path("*.out"), emit: output
  path("*_mqc.tsv"), emit: mqc
  path("versions.txt"), emit: versions

  script:
  ctrl = meta.control != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}_${meta.control}"
  """
  echo "epic2 "\$(epic2 --version) &> versions.txt
  epic2 -t ${bam} \\
    ${ctrl} \\
    --chromsizes ${chromsize} \\
    --effective-genome-fraction ${effGenomeSize} \\
    -o ${prefix}_epic.out \\
    ${args}

  echo "track type=broadPeak name=\"${meta.id}\" description=\"${prefix}\" nextItemButton=on" > ${prefix}_epic2_peaks.broadPeak
  awk -v id="${meta.id}" 'NF<10{fc=0;pval=-1;qval=-1}NF==10{fc=\$10;pval=\$4;qval=\$9}NR>1{OFS="\t"; print \$1,\$2,\$3,id"_peak_"NR,\$5,".",fc,pval,qval}' ${prefix}_epic.out >> ${prefix}_epic2_peaks.broadPeak
  cat ${prefix}_epic2_peaks.broadPeak | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${meta.id}", \$1 }' | cat $peakCountHeader - > ${prefix}_epic2_peaks.count_mqc.tsv
  """
}




/*
 * PhantomPeakQualTools QC
 */

process PPQT{
  tag "${prefix}"
  label 'ppqt'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(prefix), path(bam), path(bai)
  path sppCorrelationHeader 
  path sppNSCHeader
  path sppRSCHeader

  output:
  path '*.pdf', emit: ppqtPlot
  path '*.spp.out', emit: ppqtOutMqc
  path '*_mqc.tsv', emit: ppqtCsvMqc
  path("versions.txt")  , emit: versions

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  RUN_SPP=`which run_spp.R`
  Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="${bam}" -savp="${prefix}.spp.pdf" -savd="${prefix}.spp.Rdata" -out="${prefix}.spp.out" -p=$task.cpus
  cp $sppCorrelationHeader ${prefix}_spp_correlation_mqc.tsv
  Rscript -e "load('${prefix}.spp.Rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${prefix}_spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"
  awk -v OFS='\t' '{print "${prefix}", \$9}' ${prefix}.spp.out | cat $sppNSCHeader - > ${prefix}_spp_nsc_mqc.tsv
  awk -v OFS='\t' '{print "${prefix}", \$10}' ${prefix}.spp.out | cat $sppRSCHeader - > ${prefix}_spp_rsc_mqc.tsv
  """
}


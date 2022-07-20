/*
 * PhantomPeakQualTools QC
 */

process PPQT{
  tag "${meta.id}"
  label 'ppqt'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(bam), path(bai)
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
  Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="${bam}" -savp="${meta.id}.spp.pdf" -savd="${meta.id}.spp.Rdata" -out="${meta.id}.spp.out" -p=$task.cpus
  cp $sppCorrelationHeader ${meta.id}_spp_correlation_mqc.tsv
  Rscript -e "load('${meta.id}.spp.Rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${meta.id}_spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"
  awk -v OFS='\t' '{print "${meta.id}", \$9}' ${meta.id}.spp.out | cat $sppNSCHeader - > ${meta.id}_spp_nsc_mqc.tsv
  awk -v OFS='\t' '{print "${meta.id}", \$10}' ${meta.id}.spp.out | cat $sppRSCHeader - > ${meta.id}_spp_rsc_mqc.tsv
  """
}


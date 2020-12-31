/*
 * Peak calling & annotation QC
 */
process peakQC{
  label 'r'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/peakCalling/QC/", mode: 'copy'

  when:
  !params.skipPeakQC && params.design

  input:
  path peaks 
  path annotations 
  path peakHeader

  output:
  path "*.{txt,pdf}", emit:  macsQcOutput
  path "*.tsv"      , emit:  peakMqc

  script:
  """
  plot_macs_qc.r \\
    -i ${peaks.join(',')} \\
    -s ${peaks.join(',').replaceAll("_peaks.narrowPeak","").replaceAll("_peaks.broadPeak","")} \\
    -o ./ \\
    -p peak
  plot_homer_annotatepeaks.r \\
    -i ${annotations.join(',')} \\
    -s ${annotations.join(',').replaceAll("_annotated_peaks.txt","")} \\
    -o ./ \\
    -p annotatePeaks
  cat $peakHeader annotatePeaks.summary.txt > annotatedPeaks.summary_mqc.tsv
  """
}


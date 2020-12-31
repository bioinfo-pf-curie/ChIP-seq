/************************************
 * Peaks Annotation
 */

process peakAnnoHomer{
  tag "${sampleID}"
  label 'homer'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling/annotation/", mode: 'copy'

  when:
  !params.skipPeakAnno

  input:
  tuple val(group), val(replicate), val(peaktype), val(sampleID), path (peakfile)
  path gtfFile
  path fastaFile

  output:
  path "*.txt", emit: homerMqc

  script:
  """
  annotatePeaks.pl $peakfile \\
        $fastaFile \\
        -gtf $gtfFile \\
        -cpu ${task.cpus} \\
        > ${sampleID}_annotated_peaks.txt
  """
}



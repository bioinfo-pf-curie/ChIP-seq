/************************************
 * Peaks Annotation with HOMER
 */

process annotPeaksHomer {
  tag "${prefix}"
  label 'homer'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(prefix), path(peaks)
  path gtf
  path fasta

  output:
  tuple val(prefix), path("*.txt"), emit: output

  script:
  """
  annotatePeaks.pl ${peaks} \\
        ${fasta} \\
        -gtf ${gtf} \\
        -cpu ${task.cpus} \\
        > ${peaks.baseName}_annotHomer.txt
  """
}



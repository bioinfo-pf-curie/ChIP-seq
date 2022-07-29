/*
 * Peaks Annotation with HOMER
 */

process annotatePeaks {
  tag "${meta.id}"
  label 'homer'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(peaks)
  path gtf
  path fasta

  output:
  tuple val(meta), path("*.txt"), emit: output

  script:
  """
  annotatePeaks.pl ${peaks} \\
        ${fasta} \\
        -gtf ${gtf} \\
        -cpu ${task.cpus} \\
        > ${peaks.baseName}_annotHomer.txt
  """
}



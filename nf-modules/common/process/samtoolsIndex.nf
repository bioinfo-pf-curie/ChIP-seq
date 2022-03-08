/*
 * Samtools - Index
 */

process samtoolsIndex {
  tag "$prefix"
  label 'samtools'
  label 'minCpu'
  label 'lowMem'
 
  input:
  tuple val(prefix), path (bam)

  output:
  tuple val(prefix), path("*bam.bai"), emit: bai
  path("versions.txt") , emit: versions

  script:
  """
  echo \$(samtools --version | head -1) > versions.txt
  samtools index ${bam}
  """
}
    
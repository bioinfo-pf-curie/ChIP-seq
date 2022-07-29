/*
 * Samtools - Flagstat
 */

process samtoolsFlagstat {
  tag "${meta.id}"
  label 'samtools'
  label 'minCpu'
  label 'lowMem'
 
  input:
  tuple val(meta), path (bam)

  output:
  tuple val(meta), path("*flagstats"), emit: stats
  path("versions.txt") , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${bam.baseName}"  
  """
  echo \$(samtools --version | head -1) > versions.txt
  samtools flagstat ${bam} > ${prefix}.flagstats
  """
}
    
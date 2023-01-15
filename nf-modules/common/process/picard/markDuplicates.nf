/*
 * Picard MarkDuplicates
 */

process markDuplicates {
  tag "${meta.id}"
  label 'picard'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path('*markDups.bam'), emit: bam
  path('*markDups_metrics.txt'), emit: metrics
  path('versions.txt'), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def javaArgs = task.ext.args ?: ''
  markdupMemOption = "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
  """
  echo \$(picard MarkDuplicates --version 2>&1 | sed -e 's/Version:/picard /') > versions.txt
  picard ${markdupMemOption} ${javaArgs} MarkDuplicates \\
      MAX_RECORDS_IN_RAM=50000 \\
      INPUT=${bam} \\
      OUTPUT=${prefix}_markDups.bam \\
      METRICS_FILE=${prefix}_markDups_metrics.txt \\
      REMOVE_DUPLICATES=false \\
      ASSUME_SORTED=true \\
      PROGRAM_RECORD_ID='null' \\
      VALIDATION_STRINGENCY=LENIENT
  """
}

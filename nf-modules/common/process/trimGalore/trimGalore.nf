/*
 * TrimGalore
 */


process trimGalore {
  tag "${meta.id}"
  label 'trimgalore'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("*trimmed*fastq.gz"), emit: fastq
  tuple val(meta), path("*trimming_report.txt"), emit: logs
  path ("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''

  if (meta.singleEnd) {
  """
  [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
  trim_galore \
    $args \
    --cores ${task.cpus} \
    --basename ${prefix} \
    --gzip \
    ${prefix}.fastq.gz 2> ${prefix}_trimgalore.log
  mv ${prefix}_trimmed.fq.gz ${prefix}_trimmed.fastq.gz
  echo "trim-galore "\$(trim_galore --version 2>&1 | grep version | sed 's/^.*version //; s/Last.*\$//') > versions.txt
  """
  }else{
  """
  [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
  [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
  trim_galore \
    $args \
    --cores ${task.cpus} \
    --basename ${prefix} \
    --gzip \
    --paired \
    ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz 2> ${prefix}_trimgalore.log
  mv ${prefix}_val_1.fq.gz ${prefix}_trimmed_R1.fastq.gz
  mv ${prefix}_val_2.fq.gz ${prefix}_trimmed_R2.fastq.gz
  echo "trim-galore "\$(trim_galore --version 2>&1 | grep version | sed 's/^.*version //; s/Last.*\$//') > versions.txt 
  """
  }
}
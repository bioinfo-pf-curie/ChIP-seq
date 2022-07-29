/*
 * FastQC - Quality controls on raw reads
 */

process fastqc {
  tag "${meta.id}"
  label 'fastqc'
  label 'medCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(reads)

  output:
  path("*_fastqc.{zip,html}"), emit: results
  path("versions.txt")       , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  if (meta.singleEnd){
    """
    echo \$(fastqc --version) > versions.txt
    [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
    fastqc -q ${prefix}.fastq.gz --threads ${task.cpus} ${prefix}.fastq.gz
    """
  }else{
    """
    echo \$(fastqc --version) > versions.txt
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    fastqc -q --threads ${task.cpus} ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz
    """
  }
}

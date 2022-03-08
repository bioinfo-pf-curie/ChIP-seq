/*
 * FastQC - Quality controls on raw reads
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
 */

process fastqc {
  tag "${prefix}"
  label 'fastqc'
  label 'medCpu'
  label 'lowMem'

  input:
  tuple val(prefix), path(reads)

  output:
  path("*_fastqc.{zip,html}"), emit: results
  path("versions.txt")       , emit: versions

  script:
  if (params.singleEnd){
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

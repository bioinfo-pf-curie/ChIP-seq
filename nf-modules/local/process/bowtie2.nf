/*
 * Alignment on reference genome with Bowtie2
 * External parameters :
 * @params.singleEnd : is single-end sequencing ?
 * @params.bowtie2options : addition Bowtie2 parameters
 */

process bowtie2{
  tag "${prefix}"
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
              saveAs: {filename ->
              if (filename.indexOf(".log") > 0) "logs/$filename"
              else if (params.saveAlignedIntermediates) filename}

  input:
  //tuple val(prefix), path(reads), path(index), val(genomeName)
  tuple val(prefix), path(reads)
  path(inde)

  output:
  //tuple val(prefix), val(genomeName), path("*.bam"), emit: bam 
  tuple val(prefix), path("*.bam"), emit: bam
  path("*.log"), emit: logs
  path("versions.txt"), emit: versions

  script:
  //prefix = genomeName == genomeRef ? sample : sample + '_spike'
  inputOpts = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  """
  localIndex=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
  refName=`basename \${localIndex}`
  echo \$(bowtie2 --version | awk 'NR==1{print "bowtie2 "\$3}') > versions.txt
  bowtie2 -p ${task.cpus} \
          ${params.bowtie2Options} \
           -x \${localIndex} \
          $inputOpts > ${prefix}_\${refName}.bam 2> ${prefix}_bowtie2.log
  """
}



/*
 * Alignment on reference genome
 *
/* BOWTIE2 */
process bowtie2{
  tag "${sample} on ${genomeName}"
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
              saveAs: {filename ->
              if (filename.indexOf(".log") > 0) "logs/$filename"
              else if (params.saveAlignedIntermediates) filename}
  when:
  params.aligner == "bowtie2" && !params.inputBam

  input:
  tuple val(sample), path(reads), path(index), val(genomeBase), val(genomeName)

  output:
  tuple val(sample), path("*.bam"), emit: bam 
  path "*.log"                    , emit: mqc 
  path "v_bowtie2.txt"            , emit: version

  script:
  prefix = genomeName == params.genome ? sample : sample + '_spike'
  readCommand = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  opts = params.bowtie2Opts
  """
  bowtie2 --version &> v_bowtie2.txt
  bowtie2 -p ${task.cpus} \
          ${opts} \
           -x ${index}/${genomeBase} \
          $readCommand > ${prefix}.bam 2> ${prefix}_bowtie2.log
  """
}



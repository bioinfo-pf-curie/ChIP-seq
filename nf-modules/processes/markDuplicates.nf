/*
 * Marking duplicates
 */

process markDuplicates{
  tag "${prefix}"
  label 'picard'
  label 'lowCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/mapping", mode: 'copy',
    saveAs: {filename ->
             if (!filename.endsWith(".bam") && !filename.endsWith(".bam.bai") && params.saveAlignedIntermediates ) "stats/$filename"
             else if ( (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) && params.saveAlignedIntermediates ) filename
             else null
            }

  input:
  tuple val(prefix), path(sortedBams) 

  output:
  tuple val(prefix), path("*marked.{bam,bam.bai}"), emit: bams
  tuple val(prefix), path("*marked.flagstat")     , emit: flagstat 
  path "*marked.{idxstats,stats}"                 , emit: stats 
  path "*metrics.txt"                             , emit: picstats 
  path("v_picard.txt")                            , emit: version

  script:
  """
  echo \$(picard MarkDuplicates --version 2>&1) &> v_picard.txt
  picard -Xmx4g MarkDuplicates \\
    INPUT=${sortedBams[0]} \\
    OUTPUT=${prefix}_marked.bam \\
    ASSUME_SORTED=true \\
    REMOVE_DUPLICATES=false \\
    METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
    VALIDATION_STRINGENCY=LENIENT \\
    TMP_DIR=tmpdir
  samtools index ${prefix}_marked.bam
  samtools idxstats ${prefix}_marked.bam > ${prefix}_marked.idxstats
  samtools flagstat ${prefix}_marked.bam > ${prefix}_marked.flagstat
  samtools stats ${prefix}_marked.bam > ${prefix}_marked.stats
  """
}


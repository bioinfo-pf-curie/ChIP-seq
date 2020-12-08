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
  set val(prefix), file(sortedBams) 

  output:
  set val(prefix), file("*marked.{bam,bam.bai}")
  set val(prefix), file("*marked.flagstat") 
  file "*marked.{idxstats,stats}" 
  file "*metrics.txt" 
  file("v_picard.txt")

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


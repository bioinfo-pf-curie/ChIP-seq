/*
 * Alignment on reference genome
 */
/* STAR */
process star{
  tag "${sample} on ${genomeBase}"
  label 'star'
  label 'highCpu'
  label 'extraMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
             saveAs: {filename ->
             if (filename.indexOf(".log") > 0) "logs/$filename"
             else if (params.saveAlignedIntermediates) filename}
  when:
  params.aligner == "star" && !params.inputBam

  input:
  tuple val(sample), file(reads), file(index), val(genomeBase)

  output:
  tuple val(sample), path("*.bam"), emit: bam 
  path "*Log.final.out"           , emit: mqc 
  path "v_star.txt"               , emit: version

  script:
  prefix = genomeBase == params.genome ? sample : sample + '_spike'
  //prefix = genomeBase == genomeRef ? sample : sample + '_spike'
  opts = params.starOpts
  """
  STAR --version &> v_star.txt
  STAR --genomeDir $index \
       --readFilesIn $reads \
       --runThreadN ${task.cpus} \
       --runMode alignReads \
       --outSAMtype BAM Unsorted \
       --readFilesCommand zcat \
       --runDirPerm All_RWX \
       --outSAMunmapped Within \
       --outTmpDir /local/scratch/rnaseq_\$(date +%d%s%S%N) \
       --outFileNamePrefix $prefix  \
       --outSAMattrRGline ID:$prefix SM:$prefix LB:Illumina PL:Illumina \
       ${opts}
  """
}


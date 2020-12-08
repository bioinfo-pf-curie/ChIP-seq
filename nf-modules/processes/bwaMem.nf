/*
 * Alignment on reference genome
 */

/* BWA-MEM */
process bwaMem{
  tag "${sample} on ${genomeBase}"
  label 'bwa'
  label 'highCpu'
  label 'highMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
             saveAs: {filename ->
             if (filename.indexOf(".log") > 0) "logs/$filename"
             else if (params.saveAlignedIntermediates) filename}

  when:
  params.aligner == "bwa-mem" && !params.inputBam

  input:
  set val(sample), file(reads), file(index), val(genomeBase)

  output:
  set val(sample), file("*.bam") 
  file("*.log") 
  file("v_bwa.txt")

  script:
  prefix = genomeBase == genomeRef ? sample : sample + '_spike'
  opts = params.bwaOpts
  """
  echo \$(bwa 2>&1) &> v_bwa.txt
  bwa mem -t ${task.cpus} \
           ${index}/${genomeBase} \
          ${opts} \
          $reads | samtools view -bS - > ${prefix}.bam
  getBWAstats.sh -i ${prefix}.bam -p ${task.cpus} > ${prefix}_bwa.log
  """
}


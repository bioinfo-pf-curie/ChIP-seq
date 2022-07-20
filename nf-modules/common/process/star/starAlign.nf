/*
 * STAR reads alignment
 */

process starAlign {
  tag "$meta.id"
  label 'star'
  label 'highCpu'
  label 'extraMem'

  input:
  tuple val(meta), path(reads)
  path index
  path gtf

  output:
  tuple val(meta), path('*Aligned.out.bam'), emit: bam
  path ("*out"), emit: logs
  path ("versions.txt"), emit: versions
  tuple val(meta), path("*ReadsPerGene.out.tab"), optional: true, emit: counts
  path("*out.tab"), optional: true, emit: countsLogs
  tuple val(meta), path("*Aligned.toTranscriptome.out.bam"), optional: true, emit: transcriptsBam

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def gtfOpts = gtf.size() > 0 ? "--sjdbGTFfile ${gtf}" : ""
  """
  echo "STAR "\$(STAR --version 2>&1) > versions.txt
  STAR --genomeDir $index \\
       --readFilesIn $reads  \\
       --runThreadN ${task.cpus} \\
       --runMode alignReads \\
       --outSAMtype BAM Unsorted  \\
       --readFilesCommand zcat \\
       --runDirPerm All_RWX \\
       --outTmpDir "${params.tmpDir}/star_\$(date +%d%s%S%N)"\\
       --outFileNamePrefix $meta.id  \\
       --outSAMattrRGline ID:$meta.id SM:$meta.id LB:Illumina PL:Illumina  \\
       --outSAMunmapped Within \\
       ${gtfOpts} \\
       ${args}
  """
}

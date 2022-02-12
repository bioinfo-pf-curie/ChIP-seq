/*
 * STAR reads alignment
 */

process starAlign {
  tag "$prefix"
  label 'star'
  label 'highCpu'
  label 'extraMem'

  input:
  tuple val(prefix), path(reads)
  path index

  output:
  tuple val(prefix), path('*Aligned.out.bam'), emit: bam
  path ("*out"), emit: logs
  path ("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
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
       --outFileNamePrefix $prefix  \\
       --outSAMattrRGline ID:$prefix SM:$prefix LB:Illumina PL:Illumina  \\
       --outSAMunmapped Within \\
       ${args}
  """
}

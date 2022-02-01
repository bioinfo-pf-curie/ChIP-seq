/*
 * STAR reads alignment
 * External parameters :
 * @ params.starAlignOptions : addition STAR alignment parameters
 */

process starAlign {
  tag "$prefix"
  label 'star'
  label 'highCpu'
  label 'extraMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf(".bam") == -1) "logs/$filename"
      else if (params.saveAlignedIntermediates) filename
      else null
    }
  publishDir "${params.outDir}/counts", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("ReadsPerGene.out.tab") > 0) "$filename"
      else null
    }

  input:
  tuple val(prefix), path(reads)
  path index
  path gtf
  val useGtf

  output:
  tuple val(prefix), path('*Aligned.out.bam'), emit: bam
  tuple val(prefix), path("*ReadsPerGene.out.tab"), optional: true, emit: counts
  path("*out.tab"), optional: true, emit: countsLogs
  tuple val(prefix), path("*Aligned.toTranscriptome.out.bam"), optional: true, emit: transcriptsBam
  path ("*out"), emit: logs
  path ("versions.txt"), emit: versions

  script:
  gtfOpts = useGtf ? "--sjdbGTFfile ${gtf}" : ""
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
       ${gtfOpts} \\
       ${params.starAlignOptions}
  """
}

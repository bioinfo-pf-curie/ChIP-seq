/*
 * BAM Filtering
 */

process bamFiltering {
  tag "${prefix}"
  label 'samtools'
  label 'lowCpu'
  label 'lowMem'
  publishDir path: "${params.outDir}/mapping", mode: 'copy',
    saveAs: {filename ->
             if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "stats/$filename"
             else if (filename.endsWith("_filtered.bam") || (filename.endsWith("_filtered.bam.bai"))) filename
             else null}

  input:
  tuple val(prefix), path(markedBam)

  output:
  tuple val(prefix), path("*filtered.{bam,bam.bai}"), emit: filteredBams
  tuple val(prefix), path("*filtered.flagstat")     , emit: filteredFlagstat
  path "*filtered.{idxstats,stats}"                 , emit: stats
  path("v_samtools.txt")                            , emit: version

  script:
  filterParams = params.singleEnd ? "-F 0x004" : params.keepSingleton ? "-F 0x004" : "-F 0x004 -F 0x0008 -f 0x001"
  dupParams = params.keepDups ? "" : "-F 0x0400"
  mapqParams = params.mapq > 0 ? "-q ${params.mapq}" : ""
  nameSortBam = params.singleEnd ? "" : "samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $prefix ${prefix}_filtered.bam"
  """
  samtools --version &> v_samtools.txt
  samtools view \\
    $mapqParams \\
    -b ${markedBam[0]} > ${prefix}_filtered.bam
  samtools index ${prefix}_filtered.bam
  samtools flagstat ${prefix}_filtered.bam > ${prefix}_filtered.flagstat
  samtools idxstats ${prefix}_filtered.bam > ${prefix}_filtered.idxstats
  samtools stats ${prefix}_filtered.bam > ${prefix}_filtered.stats
  $nameSortBam
  """
} 

/*
 * Sorting BAM files
 */

process bamSort{
  tag "${prefix}"
  label 'samtools'
  label 'lowCpu'
  label 'lowMem'
  publishDir path: "${params.outDir}/mapping", mode: 'copy',
    saveAs: {filename ->
             if ( filename.endsWith("stats") && params.saveAlignedIntermediates ) "stats/$filename"
             else if ( (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) && params.saveAlignedIntermediates ) filename
             else null
            }

  input:
  set val(prefix), file(unsortedBam)

  output:
  set val(prefix), file('*sorted.{bam,bam.bai}')
  file("*stats") 
  file("*mqc") 
  file("v_samtools.txt")

  script:
  """
  samtools --version &> v_samtools.txt
  samtools sort $unsortedBam -@ ${task.cpus} -T ${prefix} -o ${prefix}_sorted.bam
  samtools index ${prefix}_sorted.bam
  samtools flagstat ${prefix}_sorted.bam > ${prefix}_sorted.flagstats
  samtools idxstats ${prefix}_sorted.bam > ${prefix}_sorted.idxstats
  samtools stats ${prefix}_sorted.bam > ${prefix}_sorted.stats

  aligned="\$(samtools view -@ $task.cpus -F 0x100 -F 0x4 -F 0x800 -c ${prefix}_sorted.bam)"
  hqbam="\$(samtools view -@ $task.cpus -F 0x100 -F 0x800 -F 0x4 -q 10 -c ${prefix}_sorted.bam)"
  lqbam="\$((\$aligned - \$hqbam))"
  echo -e "Mapped,\${aligned}" > ${prefix}_mappingstats.mqc
  echo -e "HighQual,\${hqbam}" >> ${prefix}_mappingstats.mqc
  echo -e "LowQual,\${lqbam}" >> ${prefix}_mappingstats.mqc
  """
}

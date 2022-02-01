/*
 * BigWig Tracks
 */

process bigWig {
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'lowMem'

  input:
  tuple val(prefix), path(bam), path(bai)
  path(BLbed)

  output:
  tuple val(prefix), path('*.bigwig'), emit: bigWig
  path("versions.txt"), emit: versions

  script:
  if (params.singleEnd){
    extend = params.fragmentSize > 0 && !params.noReadExtension ? "--extendReads ${params.fragmentSize}" : ""
  }else{
    extend = params.noReadExtension ? "" : "--extendReads"
  }
  blacklistParams = params.blacklist ? "--blackListFileName ${BLbed}" : ""
  effGsize = params.effGenomeSize ? "--effectiveGenomeSize ${params.effGenomeSize}" : ""
  """
  echo \$(bamCoverage --version ) > versions.txt
  bamCoverage --version &> v_deeptools.txt
  nbreads=\$(samtools view -@ $task.cpus -F 0x100 -F 0x4 -F 0x800 -c ${bam})
  sf=\$(echo "1000000 \$nbreads" | awk '{print \$1/\$2}')

  bamCoverage -b ${bam} \\
              -o ${prefix}_norm.bigwig \\
              -p ${task.cpus} \\
              ${blacklistParams} \\
              ${effGsize} \\
              ${extend} \\
              --scaleFactor \$sf
  """
}


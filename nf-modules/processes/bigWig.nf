/*
 * BigWig Tracks
 */

process bigWig {
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'lowMem'
  publishDir "${params.outDir}/bigWig", mode: "copy",
    saveAs: {filename ->
             if ( filename.endsWith(".bigwig") ) "$filename"
             else null}

  input:
  tuple val(prefix), path(filteredBams)
  path(BLbed)

  output:
  tuple val(prefix), path('*.bigwig'), emit: bigWig
  path("v_deeptools.txt")          , emit: version

  script:
  if (params.singleEnd){
    extend = params.fragmentSize > 0 && !params.noReadExtension ? "--extendReads ${params.fragmentSize}" : ""
  }else{
    extend = params.noReadExtension ? "" : "--extendReads"
  }
  blacklistParams = params.blacklist ? "--blackListFileName ${BLbed}" : ""
  effGsize = params.effGenomeSize ? "--effectiveGenomeSize ${params.effGenomeSize}" : ""
  """
  bamCoverage --version &> v_deeptools.txt
  nbreads=\$(samtools view -@ $task.cpus -F 0x100 -F 0x4 -F 0x800 -c ${filteredBams[0]})
  sf=\$(echo "10000000 \$nbreads" | awk '{printf "%.2f", \$1/\$2}')

  bamCoverage -b ${filteredBams[0]} \\
              -o ${prefix}_norm.bigwig \\
              -p ${task.cpus} \\
              ${blacklistParams} \\
              ${effGsize} \\
              ${extend} \\
              --scaleFactor \$sf
  """
}


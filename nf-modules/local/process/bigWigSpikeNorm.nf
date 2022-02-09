process bigWigSpikeNorm{
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'medMem'
        
  input:
  tuple val(prefix), path(bam), path(bai), val(normFactor)
  path(bed)
    
  output:
  tuple val(prefix), path('*.bigwig'), emit: bigwig
  path("versions.txt"), emit: versions    

  script:
  if (params.singleEnd){
    extend = params.fragmentSize > 0 && !params.noReadExtension ? "--extendReads ${params.fragmentSize}" : ""
  }else{
    extend = params.noReadExtension ? "" : "--extendReads"
  } 
  blacklistParams = params.blacklist ? "--blackListFileName ${bed}" : ""
  effGsize = params.effGenomeSize ? "--effectiveGenomeSize ${params.effGenomeSize}" : ""
  """
  echo \$(bamCoverage --version ) > versions.txt
  bamCoverage -b ${bam} \\
              -o ${prefix}_spikenorm.bigwig \\
              -p ${task.cpus} \\
               ${blacklistParams} \\
               ${effGsize} \\
               ${extend} \\
              --scaleFactor ${normFactor}
  """         
}

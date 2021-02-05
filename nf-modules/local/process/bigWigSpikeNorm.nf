process bigWigSpikeNorm{
    tag "${prefix}"
    label 'deeptools'
    label 'medCpu'
    label 'medMem'
    publishDir "${params.outDir}/bigWigSpike", mode: "copy",
      saveAs: {filename ->
        if ( filename.endsWith(".bigwig") ) "$filename"
        else null}
        
    input:
    tuple val(prefix), path(filteredBams), val(normFactor)
    path(BLbed)
    
    output:
    tuple val(prefix), path('*.bigwig'), emit: bigWigSF
    
    script:
    if (params.singleEnd){
      extend = params.fragmentSize > 0 && !params.noReadExtension ? "--extendReads ${params.fragmentSize}" : ""
    }else{
      extend = params.noReadExtension ? "" : "--extendReads"
    } 
    blacklistParams = params.blacklist ? "--blackListFileName ${BLbed}" : ""
    effGsize = params.effGenomeSize ? "--effectiveGenomeSize ${params.effGenomeSize}" : ""
    """
    bamCoverage -b ${filteredBams[0]} \\
                -o ${prefix}_spikenorm.bigwig \\
                -p ${task.cpus} \\
                 ${blacklistParams} \\
                 ${effGsize} \\
                 ${extend} \\
                --scaleFactor ${normFactor}
    """         
  }

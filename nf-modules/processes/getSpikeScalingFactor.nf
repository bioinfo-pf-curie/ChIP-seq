process getSpikeScalingFactor {
    label 'r'
    label 'minCpu'
    label 'medMem'
    publishDir "${params.outDir}/bigWigSpike", mode: "copy"

    input:
    path(tab)

    output:
    path "*.sf", emit: tabSF

    script:
    """
    getDESeqSF.r ${tab}
    """
  }

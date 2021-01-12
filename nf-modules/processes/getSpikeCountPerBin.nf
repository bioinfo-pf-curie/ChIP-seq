process getSpikeCountPerBin {
    label 'deeptools'
    label 'medCpu'
    label 'medMem'
    publishDir "${params.outDir}/bigWigSpike", mode: "copy"

    input:
    path(allBams)
    path(allBai)

    output:
    path "readCounts_spike_10kbins.tab"

    script:
    """
    multiBamSummary bins \\
                   -b $allBams \\
                   --binSize 10000 \\
                   -o results.npz \\
                   --outRawCounts readCounts_spike_10kbins.tab
    """
  }


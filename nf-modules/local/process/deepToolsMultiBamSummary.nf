/*
 * Deeptools - multiBamSummary
 */

process deepToolsMultiBamSummary {
  label 'deeptools'
  label 'medCpu'
  label 'medMem'
  
  input:
  path(bams)
  path(bai)

  output:
  path("readCounts_spike_10kbins.tab"), emit: output
  path("versions.txt"), emit: versions

  script:
  """
  echo \$(deeptools --version ) > versions.txt
  multiBamSummary bins \\
                 -b $bams \\
                 --binSize 10000 \\
                 -o results.npz \\
                 --outRawCounts readCounts_spike_10kbins.tab
  """
}


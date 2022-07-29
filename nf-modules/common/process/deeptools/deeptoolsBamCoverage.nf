/*
 * BigWig tracks from BAM file
 */

process deeptoolsBamCoverage {
  tag "${meta.id}"
  label 'deeptools'
  label 'medCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam), path(bai), val(sf)
  path(blacklistBed)
  val(effGenomeSize)

  output:
  tuple val(meta), path('*.bigwig'), emit: bigwig
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  blacklistOpts = blacklistBed.size() ? "--blackListFileName ${blacklistBed}" : ""
  effGsizeOpts = effGenomeSize.size() ? "--effectiveGenomeSize ${effGenomeSize[0]}" : ""
  sfOpts = sf != null ? "--scaleFactor $sf" : ""
  strandOpts = meta.strandness == 'forward' ? '--filterRNAstrand forward' : meta.strandness == 'reverse' ? '--filterRNAstrand reverse' : ''

  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(bamCoverage --version ) > versions.txt

  bamCoverage -b ${bam} \\
              -o ${prefix}.bigwig \\
              -p ${task.cpus} \\
              ${blacklistOpts} \\
              ${effGsizeOpts} \\
              ${args} \\
              ${sfOpts} \\
	      ${strandOpts}
  """
}
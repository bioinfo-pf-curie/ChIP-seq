/*
 * Calculate scaling factor using sequencing depth
 */

process getDepthNormFactor {
  tag "${meta.id}"
  label 'deeptools'
  label 'medCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), stdout, emit: sf
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(samtools --version | head -1 ) > versions.txt
  nbreads=\$(samtools view -@ $task.cpus -F 0x100 -F 0x4 -F 0x800 -c ${bam})
  sf=\$(echo "1000000 \$nbreads" | awk '{print \$1/\$2}')
  echo -e "\$sf"
  """
}


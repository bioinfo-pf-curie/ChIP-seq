/*
 * Alignement on reference genome with Bwa-mem
 */

process bwaMem{
  tag "${meta.id}"
  label 'bwa'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(reads)
  path(index)

  output:
  tuple val(meta), path("*.bam"), emit: bam
  path("*.log"), emit: logs
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  localIndex=`find -L ./ -name "*.amb" | sed 's/.amb//'`
  refName=`basename \${localIndex}`

  bwa \
    mem \
    $args \
    -t $task.cpus \
    \${localIndex} \
    $reads | samtools view -bS -@ $task.cpus -o ${prefix}_\${refName}.bam -

  getBWAstats.sh -i ${prefix}_\${refName}.bam -p ${task.cpus} > ${prefix}_bwa.log
  echo "Bwa-mem "\$(bwa 2>&1 | grep Version | cut -d" " -f2) &> versions.txt

  """
}

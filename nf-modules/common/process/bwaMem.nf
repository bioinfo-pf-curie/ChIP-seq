/*
 * Alignement on reference genome with Bwa-mem
 */

process bwaMem{
  tag "${prefix}"
  label 'bwa'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(prefix), path(reads)
  path(index)

  output:
  tuple val(prefix), path("*.bam"), emit: bam 
  path("*.log"), emit: logs
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  """
  localIndex=`find -L ./ -name "*.amb" | sed 's/.amb//'`
  refName=`basename \${localIndex}`
  echo "Bwa-mem "\$(bwa 2>&1 | grep Version | cut -d" " -f2) &> versions.txt
  bwa mem -t ${task.cpus} \
           \${localIndex} \
          ${args} \
          $reads | samtools view -bS - > ${prefix}_\${refName}.bam
  getBWAstats.sh -i ${prefix}_\${refName}.bam -p ${task.cpus} > ${prefix}_bwa.log
  """
}


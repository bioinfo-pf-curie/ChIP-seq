/*
 * DeepTools Compute Matrix
 */

process deeptoolsComputeMatrix{
  tag "${meta.id}"
  label 'deeptools'
  label 'medCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bigwig)
  path geneBed 

  output:
  path("*.{mat,gz,tab,pdf}"), emit: output
  path("*.plotProfile.tab"), emit: mqc
  path("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''
  """
  echo \$(deeptools --version ) > versions.txt
  computeMatrix scale-regions \\
                -R ${geneBed} \\
                -S ${bigwig} \\
                -o ${prefix}_matrix.mat.gz \\
                --outFileNameMatrix ${prefix}.computeMatrix.vals.mat \\
                ${args} \\
                -p ${task.cpus}
                
  plotProfile -m ${prefix}_matrix.mat.gz \\
              -o ${prefix}_bams_profile.pdf \\
              --outFileNameData ${prefix}.plotProfile.tab \\
              ${args2}
 
  ##sed -e 's/.0\t/\t/g' ${prefix}.plotProfile.tab | sed -e 's@.0\$@@g' > ${prefix}_plotProfile_mqc.tsv
  """
} 


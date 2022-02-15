/*
 * DeepTools QC
 */

process deepToolsComputeMatrix{
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'lowMem'

  input:
  tuple val(prefix), path(bigwig)
  path geneBed 

  output:
  path("*.{mat,gz,tab,pdf}"), emit: output
  path("*mqc.tab"), emit: mqc
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
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
              --outFileNameData ${prefix}.plotProfile.tab
              
  sed -e 's/.0\t/\t/g' ${prefix}.plotProfile.tab | sed -e 's@.0\$@@g' > ${prefix}_plotProfile_mqc.tab
  """
} 


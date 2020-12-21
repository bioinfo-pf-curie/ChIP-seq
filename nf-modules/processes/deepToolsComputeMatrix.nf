/*
 * DeepTools QC
 */

process deepToolsComputeMatrix{
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'lowMem'
  publishDir "${params.outDir}/deepTools/computeMatrix", mode: "copy"

  when:
  !params.skipDeepTools

  input:
  tuple val(prefix), path(bigwig)
  path geneBed 

  output:
  path("*.{mat,gz,pdf}"), emit: deeptoolsSingle
  path("*mqc.tab")      , emit: deeptoolsSingleMqc

  script:
  """
  computeMatrix scale-regions \\
                -R ${geneBed} \\
                -S ${bigwig} \\
                -o ${prefix}_matrix.mat.gz \\
                --outFileNameMatrix ${prefix}.computeMatrix.vals.mat \\
                --downstream 2000 --upstream 2000 --skipZeros --binSize 100\\
                -p ${task.cpus}
                
  plotProfile -m ${prefix}_matrix.mat.gz \\
              -o ${prefix}_bams_profile.pdf \\
              --outFileNameData ${prefix}.plotProfile.tab
              
  sed -e 's/.0\t/\t/g' ${prefix}.plotProfile.tab | sed -e 's@.0\$@@g' > ${prefix}_plotProfile_mqc.tab
  """
} 


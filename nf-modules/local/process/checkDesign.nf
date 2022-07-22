/*
 * Check design and sample plan files
 */

process checkDesign{
  label 'python'
  label 'minCpu'
  label 'minMem'
  publishDir "${params.summaryDir}/", mode: 'copy'

  input:
  path design
  path samplePlan

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args ?: ''
  """
  checkDesign.py -d $design -s $samplePlan ${args}
  """
}

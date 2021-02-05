/*
 * Check design and sampleplan files
 */

process checkDesign{
  label 'python'
  label 'minCpu'
  label 'minMem'
  publishDir "${params.summaryDir}/", mode: 'copy'

  when:
  params.design

  input:
  path design
  path samplePlan

  script:
  optSE = params.singleEnd ? "--singleEnd" : ""
  """
  checkDesign.py -d $design -s $samplePlan ${optSE}
  """
}

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

  script:
  optSE = meta.singleEnd ? "--singleEnd" : ""
  """
  checkDesign.py -d $design -s $samplePlan ${optSE}
  """
}

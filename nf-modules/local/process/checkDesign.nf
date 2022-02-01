/*
 * Check design and sample plan files
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
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
  optSE = params.singleEnd ? "--singleEnd" : ""
  """
  checkDesign.py -d $design -s $samplePlan ${optSE}
  """
}

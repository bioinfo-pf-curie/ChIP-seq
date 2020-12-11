/* 
 * include requires tasks 
 */
include { checkDesign } from '../processes/checkDesign'
include { fastQC } from '../processes/fastQC'

/* 
 * define the data analysis workflow 
 */
workflow qcFlow {
    // required inputs
    take:
      design
      samplePlan
      reads 
    // workflow implementation
    main:
      checkDesign(design, samplePlan)
      fastQC(reads)
}


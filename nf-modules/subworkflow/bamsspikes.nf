/* 
 * all Spikes analysis 
 */

/* 
 * include requires tasks 
 */
include { getSpikeCountPerBin } from '../processes/getSpikeCountPerBin' 
include { getSpikeScalingFactor } from '../processes/getSpikeScalingFactor'
include { bigWigSpikeNorm } from '../processes/bigWigSpikeNorm'

workflow bamsSpikesFlow {
    // required inputs
    take:
     chBamsSpikes 
     chBamsChip
     chBlacklist
    // workflow implementation
    main:

      getSpikeCountPerBin(
         chBamsSpikes.map{it[1][0]}.collect(),
         chBamsSpikes.map{it[1][1]}.collect()
      ) | getSpikeScalingFactor()

      getSpikeScalingFactor.out.tabSF 
        .splitCsv(header:false, sep:',')
        .map { row -> [row[0], row[1]]}
        .set{chScaleFactor}

      
      chBamsChip
        .combine(chScaleFactor)
        .filter{it[0] == it[2]}
        .map { it -> it[0,1,3]}
        .set{chBigWigScaleFactor}

      bigWigSpikeNorm(
         chBigWigScaleFactor,
         chBlacklist.collect()
      ) 

     emit:
     chBigWigSF = bigWigSpikeNorm.out.bigWigSF // channel: [ tuple val(prefix), path('*.bigwig') ]
}


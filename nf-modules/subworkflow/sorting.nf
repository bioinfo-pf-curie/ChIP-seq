/****************
 * Spike-in
 * Sorting BAM files
 */

/* 
 * include requires tasks 
 */
include { compareRefSpike } from '../processes/compareRefSpike' 
include { bamSort } from '../processes/bamSort'
include { checkMappingLog } from '../functions'

workflow sortingFlow {
    // required inputs
    take:
      chAlignReads
      useSpike
    // workflow implementation
    main:
      if (useSpike){

      /* Split and rebuild Channel to be sure of order between bams */
      chAlignRef = Channel.create()
      chAlignSpike = Channel.create()
      chAlignReads.choice( chAlignSpike, chAlignRef ){ it -> it[1] =~ 'spike' ? 1 : 0 }

      chCompAln = chAlignRef
       .join(chAlignSpike)

       
      // Spike-in
      // Merging, if necessary reference aligned reads and spike aligned reads
      compareRefSpike(chCompAln)

      // Filter removes all 'aligned' channels that fail the check
      chSpikeBams = compareRefSpike.out.spikeBams
      chSpikeBams
       .filter { sample, logs, bams -> checkMappingLog(logs, t="$params.spikePercentFilter") }
       .map { row -> [row[0], row[2]]}
       .set { chSpikeCheckBams }

      // concat spike and ref bams
      chRefBams = compareRefSpike.out.refBams
      chMappingSpikeMqc = compareRefSpike.out.mappingSpikeMqc
       .concat(chSpikeCheckBams)
       .set {chAllBams}
      }else{
        chAllBams = chAlignReads
        chMappingSpikeMqc = Channel.empty()
      }

      // Sorting BAM files
      bamSort(chAllBams) 

    emit:

      chMappingSpikeMqc                 // channel: [ path("*.log") ]
      sortBams = bamSort.out.sortBams // channel: [ val(prefix), path("*sorted.{bam,bam.bai}") ]
      chStatsMqc = bamSort.out.statMqc  // channel: [ tuple("*mqc") ]
      chSamtoolsVersionBamSort = bamSort.out.version // channel: [ tuple("v_samtools.txt") ]
}


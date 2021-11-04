/****************
 * Spike-in
 * Sorting BAM files
 */

/* 
 * include requires tasks 
 */
include { compareRefSpike } from '../process/compareRefSpike' 
include { bamSort } from '../process/bamSort'
include { checkMappingLog } from '../functions'

workflow sortingFlow {
    // required inputs
    take:
      chAlignReads
      useSpike
    // workflow implementation
    main:
    if (useSpike){
      result = Channel.empty()
      /* Split and rebuild Channel to be sure of order between bams */
      chAlignReads
      .branch {
        chAlignSpike: it[1] =~ 'spike'
        chAlignRef: !(it[1] =~ 'spike')
      }
      .set { result }

      // result.chAlignSpike.view { "$it is spike" }
      // result.chAlignRef.view { "$it is no spike" }
      
      chCompAln = result.chAlignSpike
       .join(result.chAlignRef)

      // Spike-in
      // Merging, if necessary reference aligned reads and spike aligned reads
      compareRefSpike(chCompAln)
      
      // chCompAln.view { "chCompAln value: $it" }

      // Filter removes all 'aligned' channels that fail the check
      chSpikeBams = compareRefSpike.out.spikeBams
      chSpikeBams
       .filter { sample, logs, bams -> checkMappingLog(logs, t="$params.spikePercentFilter") }
       .map { row -> [row[0], row[2]]}
       .set { chSpikeCheckBams }

      chMappingSpikeMqc = compareRefSpike.out.mappingSpikeMqc
      // concat spike and ref bams
      chRefBams = compareRefSpike.out.refBams
      chRefBams
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


/****************
 * Spike-in
 */

spikes_poor_alignment = []
def checkMappingLog(logs, t='1') {
  def nb_ref = 0;
  def nb_spike = 0;
  def percent_spike = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /Reads on ref\t([\d\.]+)/)) {
      nb_ref = matcher[0][1]                                                                              
    }else if ((matcher = line =~ /Reads on spike\t([\d\.]+)/)) {
      nb_spike = matcher[0][1]                                                                            
    }
  }
  logname = logs.getBaseName() - '_ref_bamcomp.mqc'
  percent_spike = nb_spike.toFloat() / (nb_spike.toFloat() + nb_ref.toFloat()) * 100
  percent_spike = percent_spike.round(5)
  if(percent_spike.toFloat() <= t.toFloat() ){
      log.info "#################### VERY POOR SPIKE ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_spike}% <<"
      spikes_poor_alignment << logname
      return false
  }else {
      log.info "          Passed alignment ($logname)   >> ${percent_spike}% <<"
      return true
  }
}

if (useSpike){

   /* Split and rebuild Channel to be sure of order between bams */
   chAlignRef = Channel.create()
   chAlignSpike = Channel.create()
   chAlignReads.choice( chAlignSpike, chAlignRef ){ it -> it[1] =~ 'spike' ? 1 : 0 }

   chCompAln = chAlignRef
      .join(chAlignSpike)

   // Merging, if necessary reference aligned reads and spike aligned reads
   process compareRefSpike{
     tag "${sample}"
     label 'compbam'
     label 'minCpu'
     label 'medMem'
     publishDir "${params.outDir}/spike", mode: 'copy',
              saveAs: {filename ->
              if (filename.indexOf(".log") > 0) "logs/$filename"
              else filename }

     input:
     set val(sample), file(unsortedBamSpike), file(unsortedBamRef) 

     output:
     set val(sample), file('*_clean.bam') 
     set val(sampleSpike), file('*.mqc'), file('*_clean_spike.bam')
     file ('*.mqc')
     file ('*.log')

     script:
     sampleSpike = sample + '_spike'
     """
     samtools sort $unsortedBamRef -n -@ ${task.cpus} -T ${sample}_ref -o ${sample}_ref_sorted.bam
     samtools sort $unsortedBamSpike -n -@ ${task.cpus} -T ${sample}_spike -o ${sample}_spike_sorted.bam
     compareAlignments.py -a ${sample}_ref_sorted.bam -b ${sample}_spike_sorted.bam -oa ${sample}_clean.bam -ob ${sample}_clean_spike.bam 2> ${sample}_compareAlignments.log
     """
   }

  // Filter removes all 'aligned' channels that fail the check
  chSpikeBams
        .filter { sample, logs, bams -> checkMappingLog(logs, t="$params.spikePercentFilter") }
        .map { row -> [row[0], row[2]]}
        .set { chSpikeCheckBams }


  // concat spike and ref bams
  chRefBams
    .concat(chSpikeCheckBams)
    .set {chAllBams}
}else{
  chAllBams = chAlignReads
  chMappingSpikeMqc = Channel.empty()
}


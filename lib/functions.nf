// Analysis Pipeline : custom functions 

/**
 * Check mapping percent from falgstat files
 */

def checkAlignmentPercent(meta, logs) {
  def percentAligned = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /([\d\.]+) \+ ([\d\.]+) mapped \s*/)) {
      nbAligned = matcher[0][1]
    } else if ((matcher = line =~ /([\d\.]+) \+ ([\d\.]+) in total \s*/)) {
      nbTotal = matcher[0][1]
    }
  }
  percentAligned = nbAligned.toFloat() / nbTotal.toFloat() * 100
  if(percentAligned.toFloat() <= '2'.toFloat() ){
      log.info "###### VERY POOR SPIKE ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($meta.id)    >> ${percentAligned}% <<"
      return false
  } else {
      log.info "          Passed alignment > ${meta.id} >> ${percentAligned}% <<"
      return true
  }
}


/**
 * Check the fraction of reads aligned on spike genome
 *
 */

def checkSpikeAlignmentPercent(meta, logs, t) {
  def nbRef = 0;
  def nbSpike = 0;
  def percentSpike = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /Reads on ref\t([\d\.]+)/)) {
      nbRef = matcher[0][1]
    }else if ((matcher = line =~ /Reads on spike\t([\d\.]+)/)) {
      nbSpike = matcher[0][1]
    }
  }
  percentSpike = nbSpike.toFloat() / (nbSpike.toFloat() + nbRef.toFloat()) * 100
  percentSpike = percentSpike.round(5)
  if( percentSpike.toFloat() <= t.toFloat() ){
    log.info("###### TOO LITTLE SPIKES DETECTED! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($meta.id) >> ${percentSpike}% <<")
    return false
  }else {
    log.info("###### Passed spike-in alignment ($meta.id) >> ${percentSpike}% <<")
    return true
  }
}


/**
 * Load ChIP-seq design file
 *
 * @params design
 */

def loadDesign(design){
  return Channel
    .fromPath(design)
    .ifEmpty { exit 1, "Design file not found: ${design}" }
    .splitCsv(header:true)
    .map { row ->
      if(row.CONTROLID==""){row.CONTROLID='NO_INPUT'}
      return [ row.SAMPLEID, row.CONTROLID, row.GROUP, row.PEAKTYPE ]
     }
}
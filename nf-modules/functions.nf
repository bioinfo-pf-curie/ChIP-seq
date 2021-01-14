/****************
 * Spike-in
 */

spikes_poor_alignment = []
def checkMappingLog(logs, t='1', spikes_poor_alignment) {
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
  spikes_poor_alignment << logname
}


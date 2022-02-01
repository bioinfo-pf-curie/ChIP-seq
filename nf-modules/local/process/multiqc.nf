/*
 * MultiQC for RNA-seq report
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
 */

process multiqc {
  label 'multiqc'
  label 'minCpu'
  label 'minMem'

  input:
  val customRunName
  path splan
  path metadata
  path multiqcConfig
  path design
  path ('fastqc/*') 
  path ('mapping/*') 
  path ('mapping/*')
  path ('mapping/*') 
  path ('filtering/*')
  path ('filtering/*')
  path ('preseq/*')
  path ('fragSize/*') 
  path ('ppqt/*') 
  path ('ppqt/*')
  path ('deepTools/*') 
  path ("deepTools/*")
  path ("deepTools/*") 
  path ('peakCalling/*') 
  path ('peakCalling/*')
  path ('peakCalling/*')
  path ('peakQC/*')
  path ('softwareVersions/*')
  path ('workflowSummary/*')
  path warnings
 

  output:
  path splan
  path "*_report.html", emit: multiqc_report
  path "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_chipseq_report" : "--filename chipseq_report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  isPE = params.singleEnd ? "" : "-p"
  designOpts= params.design ? "-d ${params.design}" : ""
  modules_list = "-m custom_content -m fastqc -m bowtie2 -m star -m preseq -m picard -m phantompeakqualtools -m deeptools -m macs2 -m homer"
  """
  stats2multiqc.sh -s ${splan} ${designOpts} -a ${params.aligner} ${isPE} > mqc.stats
  medianReadNb="\$(sort -t, -k3,3n mqc.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) printf "%.0f", (a[x-1]+a[x])/2; else printf "%.0f",a[x-1];}')"
  mqc_header.py --splan ${splan} --name "ChIP-seq" --version ${workflow.manifest.version} ${metadataOpts} --nbreads \${medianReadNb} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c multiqc-config-header.yaml -c $multiqcConfig $modules_list
  """
}


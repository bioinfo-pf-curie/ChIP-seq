process getSoftwareVersions{
  label 'python'
  label 'minCpu'
  label 'minMem'
  publishDir path: "${params.outDir}/softwareVersions", mode: "copy"

  when:
  !params.skipSoftVersions

  input:
  path 'v_fastqc.txt' 
  path 'v_bwa.txt' 
  path 'v_bowtie2.txt' 
  path 'v_star.txt'
  path 'v_samtools.txt' 
  path 'v_picard.txt' 
  path 'v_macs2.txt'
  path 'v_epic2.txt' 
  path 'v_preseq.txt'
  path 'v_idr.txt'
  path 'v_R.txt' 
  path 'v_deeptools.txt'
  path 'v_featurecounts.txt'

  output:
  path 'software_versions_mqc.yaml'
  
  script:
  """
  echo $workflow.manifest.version &> v_pipeline.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
} 


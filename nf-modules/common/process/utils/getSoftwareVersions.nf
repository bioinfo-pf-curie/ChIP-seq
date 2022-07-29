/*
 * getSoftwareVersions
 */


process getSoftwareVersions{
  label 'python'
  label 'minCpu'
  label 'lowMem'

  input:
  path versions

  output:
  path 'software_versions_mqc.yaml', emit: versionsYaml

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  echo "Pipeline $workflow.manifest.version" > all_versions.txt
  echo "Nextflow $workflow.nextflow.version" >> all_versions.txt
  cat ${versions} >> all_versions.txt
  scrape_software_versions.py -i all_versions.txt > software_versions_mqc.yaml
  """
}

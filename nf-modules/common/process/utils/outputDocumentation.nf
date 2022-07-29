/*
 * Output documentation in HTML format
 */

process outputDocumentation {
  label 'python'
  label 'minCpu'
  label 'minMem'

  input:
  path outputDocs 
  path images 

  output:
  path "results_description.html"

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  markdown_to_html.py $outputDocs -o results_description.html
  """
}

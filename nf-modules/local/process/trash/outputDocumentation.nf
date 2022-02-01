process outputDocumentation {
    label 'python'
    label 'minCpu'
    label 'minMem'
    publishDir "${params.summaryDir}/", mode: 'copy'

    input:
    path output_docs 
    path images 

    output:
    path "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

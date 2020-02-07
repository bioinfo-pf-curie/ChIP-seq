#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    'nf-core/chipseq': ['soft_versions.txt', r"(\S+)"],
    'Nextflow': ['soft_versions.txt', r"(\S+)"],
    'FastQC': ['soft_versions.txt', r"FastQC v(\S+)"],
    'MultiQC': ['soft_versions.txt', r"multiqc, version (\S+)"],
	'BWA': ['soft_versions.txt', r"Version: (\S+)"],
	'Bowtie2': ['soft_versions.txt', r"version (\S+)"],
	'STAR': ['soft_versions.txt', r"STAR(\S+)"],
	'samtools': ['soft_versions.txt', r"samtools (\S+)"],
	'picard': ['soft_versions.txt', r"(\S+)"],
	'preseq': ['soft_versions.txt', r"Version: (\S+)"]
}
results = OrderedDict()
results['nf-core/chipseq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['BWA'] = '<span style="color:#999999;\">N/A</span>'
results['Bowtie2'] = '<span style="color:#999999;\">N/A</span>'
results['STAR'] = '<span style="color:#999999;\">N/A</span>'
results['samtools'] = '<span style="color:#999999;\">N/A</span>'
results['picard'] = '<span style="color:#999999;\">N/A</span>'
results['preseq'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'nf-core/chipseq-software-versions'
section_name: 'nf-core/chipseq Software Versions'
section_href: 'https://github.com/nf-core/chipseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")

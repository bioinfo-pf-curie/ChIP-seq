#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    'nf-core/chipseq': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
	'BWA': ['v_bwa.txt', r"Version: (\S+)"],
	'Bowtie2': ['v_bowtie2.txt', r"version (\S+)"],
	'STAR': ['v_star.txt', r"STAR(\S+)"],
	'samtools': ['v_samtools.txt', r"samtools (\S+)"],
	'bedtools': ['v_bedtools.txt', r"bedtools v(\S+)"],
	'picard': ['v_picard.txt', r"([\d\.]+)-SNAPSHOT"],
	'preseq': ['v_preseq.txt', r"Version: (\S+)"],
	'deeptools': ['v_deeptools.txt', r"plotFingerprint (\S+)"],
	'R': ['v_R.txt', r"R version (\S+)"],
	'MACS2': ['v_macs2.txt', r"macs2 (\S+)"],
	'epic2': ['v_epic2.txt', r"(\S+)"],
	'idr': ['v_idr.txt', r"IDR (\S+)"]
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
results['bedtools'] = '<span style="color:#999999;\">N/A</span>'
results['picard'] = '<span style="color:#999999;\">N/A</span>'
results['preseq'] = '<span style="color:#999999;\">N/A</span>'
results['deeptools'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'
results['MACS2'] = '<span style="color:#999999;\">N/A</span>'
results['epic2'] = '<span style="color:#999999;\">N/A</span>'
results['idr'] = '<span style="color:#999999;\">N/A</span>'



# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'Chip-seq pipeline software versions'
section_name: 'Software Versions'
section_href: 'https://gitlab.curie.fr/data-analysis/chip-seq/t'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")

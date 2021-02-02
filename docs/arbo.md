```shell
chip-seq
├── assets
│   └── -
├── bin
│   ├── checkDesign.py
│   ├── compareAlignments.py
│   ├── getBWAstats.sh
│   ├── getDESeqSF.r
│   ├── getNormFactor.py
│   ├── markdown_to_html.py
│   ├── mqc_header.py
│   ├── plot_homer_annotatepeaks.r
│   ├── plot_macs_qc.r
│   ├── __pycache__
│   │   └── scrape_software_versions.cpython-36.pyc
│   ├── replicate_idr.py
│   ├── scrape_software_versions.py
│   └── stats2multiqc.sh
├── CHANGELOG.md
├── conf
│   ├── base.config
│   ├── cluster.config
│   ├── conda.config
│   ├── docker.config
│   ├── geniac.config
│   ├── genomes.config
│   ├── multiconda.config
│   ├── multipath.config
│   ├── path.config
│   ├── process.config
│   ├── singularity.config
│   └── test.config
├── docs
│   └── -
├── environment.yml
├── geniac
├── LICENSE
├── main.nf
├── nextflow.config
├── nf-modules
│   ├── functions.nf
│   ├── processes
│   │   ├── bamFiltering.nf
│   │   ├── bamSort.nf
│   │   ├── bigWig.nf
│   │   ├── bigWigSpikeNorm.nf
│   │   ├── bowtie2.nf
│   │   ├── broadMACS2.nf
│   │   ├── bwaMem.nf
│   │   ├── checkDesign.nf
│   │   ├── compareRefSpike.nf
│   │   ├── deepToolsComputeMatrix.nf
│   │   ├── deepToolsCorrelationQC.nf
│   │   ├── deepToolsFingerprint.nf
│   │   ├── fastQC.nf
│   │   ├── featureCounts.nf
│   │   ├── getSoftwareVersions.nf
│   │   ├── getSpikeCountPerBin.nf
│   │   ├── getSpikeScalingFactor.nf
│   │   ├── IDR.nf
│   │   ├── markDuplicates.nf
│   │   ├── multiqc.nf
│   │   ├── outputDocumentation.nf
│   │   ├── peakAnnoHomer.nf
│   │   ├── peakQC.nf
│   │   ├── PPQT.nf
│   │   ├── prepareAnnotation.nf
│   │   ├── preseq.nf
│   │   ├── sharpMACS2.nf
│   │   ├── star.nf
│   │   ├── veryBroadEpic2.nf
│   │   └── workflowSummaryMqc.nf
│   └── subworkflow
│       ├── bamschip.nf
│       ├── bamsspikes.nf
│       ├── mapping.nf
│       ├── markdup.nf
│       ├── peakcalling.nf
│       ├── qc.nf
│       └── sorting.nf
├── README.md
├── recipes
│   └── -
└── test
    └── -
```

```shell
main.nf/
└── workflow:main
    └── nf-modules
        ├── local
        │   ├── process
        │   │   ├── featureCounts
        │   │   ├── getSoftwareVersions
        │   │   ├── multiqc
        │   │   ├── outputDocumentation
        │   │   ├── prepareAnnotation
        │   │   └── workflowSummaryMqc
        │   └── subworkflow
        │       ├── bamsChipFlow
        │       │   └── process
        │       │       ├── bigWig
        │       │       ├── deepToolsComputeMatrix
        │       │       ├── deepToolsCorrelationQC
        │       │       ├── deepToolsFingerprint
        │       │       └── PPQT
        │       ├── bamsSpikesFlow
        │       │   └── process
        │       │       ├── bigWigSpikeNorm
        │       │       ├── getSpikeCountPerBin
        │       │       └── getSpikeScalingFactor
        │       ├── mappingFlow
        │       │   └── process
        │       │       ├── bowtie2
        │       │       ├── bwaMem
        │       │       └── star
        │       ├── markdupFlow
        │       │   └── process
        │       │       ├── bamFiltering
        │       │       ├── markDuplicates
        │       │       └── preseq
        │       ├── peakCallingFlow
        │       │   └── process
        │       │       ├── broadMACS2
        │       │       ├── IDR
        │       │       ├── peakAnnoHomer
        │       │       ├── peakQC
        │       │       ├── sharpMACS2
        │       │       └── veryBroadEpic2
        │       ├── qcFlow
        │       │   └── process
        │       │       ├── checkDesign
        │       │       └── fastQC
        │       └── sortingFlow
        │           └── process
        │               ├── bamSort
        │               └── compareRefSpike
        └── workflow.onComplete

```

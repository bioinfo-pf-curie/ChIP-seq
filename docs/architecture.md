```shell
main.nf/
└── workflow:main
    └── nf-modules
        ├── processes
        │   ├── featureCounts
        │   ├── getSoftwareVersions
        │   ├── multiqc
        │   ├── outputDocumentation
        │   ├── prepareAnnotation
        │   └── workflowSummaryMqc
        ├── subworkflow
        │   ├── bamsChipFlow
        │   │   └── processes
        │   │       ├── bigWig
        │   │       ├── deepToolsComputeMatrix
        │   │       ├── deepToolsCorrelationQC
        │   │       ├── deepToolsFingerprint
        │   │       └── PPQT
        │   ├── bamsSpikesFlow
        │   │   └── processes
        │   │       ├── bigWigSpikeNorm
        │   │       ├── getSpikeCountPerBin
        │   │       └── getSpikeScalingFactor
        │   ├── mappingFlow
        │   │   └── processes
        │   │       ├── bowtie2
        │   │       ├── bwaMem
        │   │       └── star
        │   ├── markdupFlow
        │   │   └── processes
        │   │       ├── bamFiltering
        │   │       ├── markDuplicates
        │   │       └── preseq
        │   ├── peakCallingFlow
        │   │   └── processes
        │   │       ├── broadMACS2
        │   │       ├── IDR
        │   │       ├── peakAnnoHomer
        │   │       ├── peakQC
        │   │       ├── sharpMACS2
        │   │       └── veryBroadEpic2
        │   ├── qcFlow
        │   │   └── processes
        │   │       ├── checkDesign
        │   │       └── fastQC
        │   └── sortingFlow
        │       └── processes
        │           ├── bamSort
        │           └── compareRefSpike
        └── workflow.onComplete
```

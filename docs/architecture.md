```shell
main.nf
└── workflow:main
    ├── bamsChipFlow
    │   ├── bigWig
    │   ├── deepToolsComputeMatrix
    │   ├── deepToolsCorrelationQC
    │   ├── deepToolsFingerprint
    │   └── PPQT
    ├── bamsSpikesFlow
    │   ├── bigWigSpikeNorm
    │   ├── getSpikeCountPerBin
    │   └── getSpikeScalingFactor
    ├── featureCounts
    ├── getSoftwareVersions
    ├── mappingFlow
    │   ├── bowtie2
    │   ├── bwaMem
    │   └── star
    ├── markdupFlow
    │   ├── bamFiltering
    │   ├── markDuplicates
    │   └── preseq
    ├── multiqc
    ├── outputDocumentation
    ├── peakCallingFlow
    │   ├── broadMACS2
    │   ├── IDR
    │   ├── peakAnnoHomer
    │   ├── peakQC
    │   ├── sharpMACS2
    │   └── veryBroadEpic2
    ├── prepareAnnotation
    ├── qcFlow
    │   ├── checkDesign
    │   └── fastQC
    ├── sortingFlow
    │   ├── bamSort
    │   └── compareRefSpike
    ├── workflow.onComplete
    └── workflowSummaryMqc
```

#!/bin/bash 

echo "nextflow run /bioinfo/users/nservant/GitLab/chip-seq/main.nf \
    --samplePlan '/bioinfo/users/nservant/GitLab/chip-seq/test/samplePlan_fulldataset.csv' --singleEnd \
    --design '/bioinfo/users/nservant/GitLab/chip-seq/test/design_fulldataset.csv' \
    --spike 'dmelr6.22' \
    --genome 'hg38' --genomeAnnotationPath '/data/annotations/pipelines' \
    -profile cluster,path --globalPath '/data/kdi_prod/.kdi/project_workspace_0/982/acl/01.00/conda/chipseq/bin' \
    --outDir '/data/tmp/nservant/chipseq_fulltest/' \
    -w '/data/tmp/chipseq_fulltest/work' -resume" | qsub -N "fulltest" -l "mem=1gb,nodes=1:ppn=1"

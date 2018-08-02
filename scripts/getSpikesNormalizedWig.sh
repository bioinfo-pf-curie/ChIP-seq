#!/bin/bash

## Method from https://github.com/deeptools/deepTools/issues/509
odir=/data/kdi_prod/.kdi/project_workspace_0/1325/acl/01.00/results/spike_analysis

## List of samples
indir=/data/kdi_prod/.kdi/project_workspace_0/1325/acl/01.00/results/processing_with_spikes
SAMPLE_LIST_DROSO="${indir}/A949C07/mapping/A949C07.R1_dmel-all-chromosome-r6.21.bam ${indir}/A949C08/mapping/A949C08.R1_dmel-all-chromosome-r6.21.bam \
${indir}/A949C09/mapping/A949C09.R1_dmel-all-chromosome-r6.21.bam ${indir}/A949C04/mapping/A949C04.R1_dmel-all-chromosome-r6.21.bam \
${indir}/A949C10/mapping/A949C10.R1_dmel-all-chromosome-r6.21.bam ${indir}/A949C05/mapping/A949C05.R1_dmel-all-chromosome-r6.21.bam \
${indir}/A949C11/mapping/A949C11.R1_dmel-all-chromosome-r6.21.bam ${indir}/A949C06/mapping/A949C06.R1_dmel-all-chromosome-r6.21.bam \
${indir}/A949C12/mapping/A949C12.R1_dmel-all-chromosome-r6.21.bam"

SAMPLE_LIST_REF="${indir}/A949C07/mapping/A949C07.R1_dmel-all-chromosome-r6.21_R1_hg38.bam ${indir}/A949C08/mapping/A949C08.R1_dmel-all-chromosome-r6.21_R1_hg38.bam \
${indir}/A949C09/mapping/A949C09.R1_dmel-all-chromosome-r6.21_R1_hg38.bam ${indir}/A949C04/mapping/A949C04.R1_dmel-all-chromosome-r6.21_R1_hg38.bam \
${indir}/A949C10/mapping/A949C10.R1_dmel-all-chromosome-r6.21_R1_hg38.bam ${indir}/A949C05/mapping/A949C05.R1_dmel-all-chromosome-r6.21_R1_hg38.bam \
${indir}/A949C11/mapping/A949C11.R1_dmel-all-chromosome-r6.21_R1_hg38.bam ${indir}/A949C06/mapping/A949C06.R1_dmel-all-chromosome-r6.21_R1_hg38.bam \
${indir}/A949C12/mapping/A949C12.R1_dmel-all-chromosome-r6.21_R1_hg38.bam"


## Make count table per bins
echo -e "Generate counts table per bin ..."
/bioinfo/local/build/Centos/python/python-2.7.11/bin/multiBamSummary bins --binSize 10000 -o ${odir}/results.npz --outRawCounts ${odir}/readCounts_10kbins.tab \
    --bamfiles $SAMPLE_LIST_DROSO

## Calculate scaling factor
echo -e "Calculate DESeq scaling factor ..."
/bioinfo/local/build/Centos/R/R-3.4.0/bin/R CMD BATCH --no-save --no-restore "--args count.table='${odir}/readCounts_10kbins.tab'" getDESeqSF.R

for bam in ${SAMPLE_LIST_REF}
do
prefix=$(basename $bam | sed -e 's/_R1_hg38.bam//')
sf=$(grep $prefix ${odir}/readCounts_10kbins.sf | awk '{print $2}')
echo -e "Generate Wig file for ${prefix} [sf=$sf] ..."
/bioinfo/local/build/Centos/python/python-2.7.11/bin/bamCoverage -b ${bam} -o ${odir}/${prefix}_spikenorm.bw --scaleFactor ${sf} 2>> bamCoverage.log
done



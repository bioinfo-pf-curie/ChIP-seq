***********************************
CHANGES IN VERSION 1.0.3

NEW FEATURES

   o Add --spikePercentFilter to fix the minimum percentage of reads which should align on spike-in genome (#23)
   o Add --keepSingleton option

BUG FIXES

   o Fix bugs in sample name (#39)
   o Fix bug in spike-in reads cleaning
   o Update BWA stats for PE data. Report pairs and not reads for mqc
   o Dups % is calculated over the number of mapped reads (not total reads)

***********************************
CHANGES IN VERSION 1.0.2

SIGNIFICANT USER-VISIBLE CHANGES

  o Report high/low quality reads instead of unique/multi reads
  o Update doc
  o Use paired and single reads in duplicates level calculation

BUG FIXES

  o export plots in output documentation (#35)
  o Add dmelr6.22 annotation genome
  o Fix bug in mapping stat for bowtie2 with paired-end data

***********************************
CHANGES IN VERSION 1.0.1

SIGNIFICANT USER-VISIBLE CHANGES

  o Update getSoftwareVersion()
  o Update process.config to use standard labels
  o update .libPaths() in all R scripts
  o Fix issue #19 with --spikeFasta and --spikeBt2Index
  o Add --fragmentSize and extension parameters for deeptools functions

BUG FIXES

  o Fix bug in DESeq2 scaling factor calculation
  

***********************************
CHANGES IN VERSION 1.0.0

NEW FEATURES

  o First stable of the chip-seq pipeline



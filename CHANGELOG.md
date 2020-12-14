version-1.0.4

NEW FEATURES

   - update STAR to 2.7.6a

***********************************
version-1.0.3

NEW FEATURES

   - Add --spikePercentFilter to fix the minimum percentage of reads which should align on spike-in genome (#23)
   - Add --keepSingleton option
   - Improve checkDesign function (#22)
   - Add --noReadExtension parameters
   - Add parameters for all mapper options

BUG FIXES

   - Fix bug for Epic2 without inputs (#47)
   - Fix bugs in sample name (#39)
   - Fix bug in spike-in reads cleaning
   - Update BWA stats for PE data. Report pairs and not reads for mqc
   - Dups % is calculated over the number of mapped reads (not total reads)

***********************************
version-1.0.2

SIGNIFICANT USER-VISIBLE CHANGES

  o Report high/low quality reads instead of unique/multi reads
  o Update doc
  o Use paired and single reads in duplicates level calculation

BUG FIXES

  - export plots in output documentation (#35)
  - Add dmelr6.22 annotation genome
  - Fix bug in mapping stat for bowtie2 with paired-end data

***********************************
version-1.0.1

SIGNIFICANT USER-VISIBLE CHANGES

  - Update getSoftwareVersion()
  - Update process.config to use standard labels
  - update .libPaths() in all R scripts
  - Fix issue #19 with --spikeFasta and --spikeBt2Index
  - Add --fragmentSize and extension parameters for deeptools functions

BUG FIXES

  - Fix bug in DESeq2 scaling factor calculation
  

***********************************
version-1.0.0

NEW FEATURES

  - First stable of the chip-seq pipeline



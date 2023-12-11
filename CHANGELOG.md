***********************************
version-2.0.0

NEW FEATURES
  - DSL2 version of the ChIP-seq pipeline
  - Update conda env
  - Add `--mac2Opts` and `--epic2Opts` to change on-the-fly the peak calling options (#59)
  - Support BAM files as inputs (#13)
  - Add `--trimming` option to remove adapter sequences (#42)

SIGNIFICANT USER-VISIBLE CHANGES
  - The 'design' file template has been updated
  - Calculate insert size with Picard for PE data (#49)
  - Add spike-in % in general table (hidden by default)(#57)

***********************************
version-1.0.6

NEW FEATURES
   - Add Mouse mm39 annotation
   - Add Human hg19/hg38 base annotation
	  
*************************************
version-1.0.5
	  
NEW FEATURES
   - Add Bombyx_v4 genome annotation
   - Add dmel6.28 genome annotation
			
SIGNIFICANT USER-VISIBLE CHANGES
	
   - Use picard for insert size distribution for paired-end data (#49)
   - Calculate SF over 1M reads instead of 10M for bigwig
				  
BUG FIXES
				  
   - Bug in mapping Channel for Bombyx (#53)
   - Effective size of Bombyx updated (#53)
						
************************************
version-1.0.4
						
NEW FEATURES
						
   - update STAR to 2.7.6a
   - add S. mikatae annotation for spikes
							  
BUG FIXES
							  
   - Fix bug when there are no blacklist regions
								 
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
  - Report high/low quality reads instead of unique/multi reads
  - Update doc
  - Use paired and single reads in duplicates level calculation

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



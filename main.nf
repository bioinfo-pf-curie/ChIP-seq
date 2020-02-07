#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. 
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
                         Chip-seq
========================================================================================
 Chip-seq Analysis Pipeline.
 #### Homepage / Documentation
 https://gitlab.curie.fr/chipseq
----------------------------------------------------------------------------------------
*/

// TODO - replace all Chip-seq with the name of your pipeline

def helpMessage() {
    // TODO: Add to this help message with new command line parameters

    if ("${workflow.manifest.version}" =~ /dev/ ){
       dev_mess = file("$baseDir/assets/dev_message.txt")
       log.info dev_mess.text
    }

    log.info"""
    
    Chip-seq v${workflow.manifest.version}
    ======================================================================

    Usage:

    nextflow run Chip-seq --reads '*_R{1,2}.fastq.gz' --genome 'hg19' -profile curie
    nextflow run Chip-seq --samplePlan sample_plan --genome hg19 -profile curie

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --samplePlan                  Path to sample plan file if '--reads' is not specified
      --genome                      Name of iGenomes reference
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    Options:
      --singleEnd                   Specifies that the input is single end reads
    
	References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help messsage
if (params.help){
    helpMessage()
    exit 0
}

// TODO - Add any reference files that are needed - see igenome.conf
// Configurable reference genomes
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if ( params.fasta ){
    ch_fasta = file(params.fasta, checkIfExists: true)
	lastPath = params.fasta.lastIndexOf(File.separator)
	bwa_base = params.fasta.substring(lastPath+1)
}
else{
	exit 1, "Fasta file not found: ${params.fasta}"
}

if (params.aligner == "bwa-mem"){
	params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
	if (params.bwa_index){
		lastPath = params.fasta.lastIndexOf(File.separator)
		bwa_dir =  params.bwa_index.substring(0,lastPath+1)
		bwa_base = params.bwa_index.substring(lastPath+1)
		Channel
			.fromPath(bwa_dir, checkIfExists: true)
			.set { ch_bwa_index }
	}
	else{
		exit 1, "BWA index file not found: ${params.bwa_index}"
	}
}

if (params.aligner == "bowtie2"){
	params.bt2_index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
	if (params.bt2_index){
		lastPath = params.fasta.lastIndexOf(File.separator)
		bt2_dir =  params.bt2_index.substring(0,lastPath+1)
		bt2_base = params.bt2_index.substring(lastPath+1)
		Channel
			.fromPath(bt2_dir, checkIfExists: true)
			.set { ch_bt2_index }	}
	else{
		exit 1, "Bowtie2 index file not found: ${params.bt2_index}"
	}
}

if (params.aligner == "star"){
	params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
	if (params.star_index){
		lastPath = params.fasta.lastIndexOf(File.separator)
		star_dir =  params.star_index.substring(0,lastPath+1)
		star_base = params.star_index.substring(lastPath+1)
		Channel
			.fromPath(star_dir, checkIfExists: true)
			.set { ch_star_index }	}
	else{
		exit 1, "STAR index file not found: ${params.star_index}"
	}
}

// Other inputs 
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
if (params.gtf) { 
	ch_gtf = file(params.gtf, checkIfExists: true) 
} 
else { 
	exit 1, "GTF annotation file not specified!" 
}

params.gene_bed = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
if (params.gene_bed)  { 
	ch_gene_bed = file(params.gene_bed, checkIfExists: true) 
}

params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false
if (params.blacklist) { 
	ch_blacklist = Channel.fromPath(params.blacklist, checkIfExists: true) 
} 
else { 
	ch_blacklist = Channel.empty() 
}

// TODO - Tools option configuration - see tools.conf
// Add here the list of options that can change from a reference genome to another
// params.star_options = params.genomes[ params.genome ].star_opts ?: params.star_opts


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Stage config files
ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)
ch_output_docs = file(params.output_doc, checkIfExists: true)
if (params.singleEnd) {
    ch_bamtools_filter_config = file(params.bamtools_filter_se_config, checkIfExists: true)
} else {
    ch_bamtools_filter_config = file(params.bamtools_filter_pe_config, checkIfExists: true)
}

/*
 * CHANNELS
 */

if ( params.metadata ){
   	Channel
       	.fromPath( params.metadata )
       	.ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
       	.set { ch_metadata }
}

/*
 * Create a channel for input read files
 */

if(params.samplePlan){
	if(params.singleEnd){
		Channel
        	.from(file("${params.samplePlan}"))
        	.splitCsv(header: false)
        	.map{ row -> [ row[0], [file(row[2])]] }
        	.into { raw_reads_fastqc;raw_reads_bwa;raw_reads_bt2;raw_reads_star }
   	} 
   	else{
      	Channel
	       	.from(file("${params.samplePlan}"))
         	.splitCsv(header: false)
         	.map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
         	.into { raw_reads_fastqc;raw_reads_bwa;raw_reads_bt2;raw_reads_star }
	}
   	params.reads=false
}
else if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc;raw_reads_bwa;raw_reads_bt2;raw_reads_star } 
    } 
	else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc;raw_reads_bwa;raw_reads_bt2;raw_reads_star }
    }
} 
else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc;raw_reads_bwa;raw_reads_bt2;raw_reads_star }
}

/*
 * Make sample plan if not available
 */

if (params.samplePlan){
  ch_splan = Channel.fromPath(params.samplePlan)
}
else{
  	if (params.singleEnd){
    	Channel
       		.from(params.readPaths)
       		.collectFile() {
         	item -> ["sparams.genome ? ample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
        }
       	.set{ ch_splan }
  	}

  	else{
     	Channel
       		.from(params.readPaths)
       		.collectFile() {
         	item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       	.set{ ch_splan }
  }
}

/*
 * Make BWA index if not available
 */
if (!params.bwa_index && params.aligner == "bwa-mem"){
	process build_BWAindex{
		tag "$fasta"
		publishDir "${params.outdir}/bwa_index"

		input:
			file fasta from ch_fasta

		output:
			file "bwa_index" into ch_bwa_index

		script:
			"""
			bwa index $fasta
			mkdir bwa_index && mv ${fasta}* bwa_index
			"""
	}
}

/*
 * Make Bowtie2 index if not available
 */
if (!params.bt2_index && params.aligner == "bowtie2"){
	process build_BT2index{
		tag "$fasta"bt2_index
		publishDir "${params.outdir}/bt2_index"

		input:
			file fasta from ch_fasta

		output:
			file "bt2_index" into ch_bt2_index

		script:
			"""
			bowtie2-build $fasta ${params.genome}
			mkdir bt2_index && mv ${params.genome}* bt2_index
			"""
	}
}

/*
 * Make STAR index if not available
 */
if (!params.star_index && params.aligner == "star"){
	process build_STARindex{
		tag "$fasta"
		publishDir "${params.outdir}/star_index"

		input:
			file fasta from ch_fasta

		output:
			file "star_index" into ch_star_index

		script:
			"""
			STAR --runMode genomeGenerate --genomeFastaFiles $fasta
			mkdir star_index && mv ${fasta}* star_index
			"""
	}
}

/*
 * Make BED file for blacklisted regions
 */

if (!params.gene_bed) {
    process makeGeneBED {
        tag "$gtf"
        publishDir "${params.outdir}/reference_genome", mode: 'copy'

        input:
        file gtf from ch_gtf

        output:
        file "*.bed" into ch_gene_bed

        script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}

/*
 * Make genome filter
 */

if(!params.skip_filtering){
	process makeGenomeFilter {
		publishDir "${params.outdir}/reference_genome", mode: 'copy'

		input:
		file fasta from ch_fasta
		file blacklist from ch_blacklist.ifEmpty([])

		output:
		file "$fasta" into ch_genome_fasta                 // FASTA FILE FOR IGV
		file "*.fai" into ch_genome_fai                    // FAI INDEX FOR REFERENCE GENOME
		file "*.bed" into ch_genome_filter_regions         // BED FILE WITHOUT BLACKLIST REGIONS
		file "*.sizes" into ch_genome_sizes_bigwig         // CHROMOSOME SIZES FILE FOR BEDTOOLS

		script:
		blacklist_filter = params.blacklist ? "sortBed -i $blacklist -g ${fasta}.sizes | complementBed -i stdin -g ${fasta}.sizes" : "awk '{print \$1, '0' , \$2}' OFS='\t' ${fasta}.sizes"
		"""
		samtools faidx $fasta
		cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
		$blacklist_filter > ${fasta}.include_regions.bed
		"""
	}
}

// Header log info
if ("${workflow.manifest.version}" =~ /dev/ ){
   dev_mess = file("$baseDir/assets/dev_message.txt")
   log.info dev_mess.text
}

log.info """=======================================================

Chip-seq v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'Chip-seq'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
// TODO : Report custom parameters here
if (params.samplePlan) {
   	summary['SamplePlan']   = params.samplePlan
}
else{
   	summary['Reads']        = params.reads
}
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile

if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
 * FastQC
 */

process fastQC{
	tag "${prefix}"
	publishDir "${params.outdir}/fastQC"

	when:
		!params.skip_fastqc

	input:
		set val(prefix), file(reads) from raw_reads_fastqc

	output:
		file "*_fastqc.{zip,html}" into ch_fastqc_results

	script:
		toqc = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
		"""
		fastqc -q $reads
		mv ${toqc}_fastqc.html ${prefix}_fastqc.html
		mv ${toqc}_fastqc.zip ${prefix}_fastqc.zip
		"""
}

/*
 * Alignment
 */

// BWA-MEM
if (!params.skip_alignment && params.aligner == "bwa-mem"){
	process BWA{
		tag "${prefix}"
		publishDir "${params.outdir}/bwa_alignment"
		
		input:
			set val(prefix), file(reads) from raw_reads_bwa
			file bwa_index from ch_bwa_index.collect()

		output:
			set val(prefix), file("*.bam") into ch_aligned_reads

		script:
			"""
			bwa mem -t 2 $bwa_index/${bwa_base} $reads \\
			| samtools view -b -h -o ${prefix}.bam
			"""
	}
}

// BOWTIE2
if (!params.skip_alignment && params.aligner == "bowtie2"){
	process bowtie2{
		tag "${prefix}"
		publishDir "${params.outdir}/bowtie2_alignment"

		input:
			set val(prefix), file(reads) from raw_reads_bt2
			file bt2_index from ch_bt2_index.collect()

		output:
			set val(prefix), file("*.bam") into ch_aligned_reads

		script:
		"""
		bowtie2 -p 2 -x $bt2_index/${bt2_base} -1 ${reads[0]} -2 ${reads[1]} \\
		| samtools view -b -h -o ${prefix}.bam
		"""
	}
}

// STAR
if (!params.skip_alignment && params.aligner == "star"){
	process star{
		tag "${prefix}"
		publishDir "${params.outdir}/star_alignment"

		input:
		input:
			set val(prefix), file(reads) from raw_reads_star
			file star_index from ch_star_index.collect()


		output:
			set val(prefix), file('*.bam') into ch_aligned_reads

		script:
		"""
		STAR --genomeDir $star_index/${star_base} --runThreadN 2 --readFilesIn $reads --outSAMtype BAM Unsorted --readFilesCommand zcat
		"""
	}
}

/*
 * Sorting BAM files
 */

if (!params.skip_alignment){
	process bamSort{
		tag "${prefix}"
		publishDir path: "${params.outdir}/sorted_bams", mode: 'copy',
					saveAs: {filename ->
								if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "samtools_stats/$filename"
								else if (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) filename
								else null 
							}
						  		
		input:
			set val(prefix), file(unsorted_bam) from ch_aligned_reads

		output:
			set val(prefix), file('*.sorted.{bam,bam.bai}') into ch_sort_bam
			file "*.{flagstat,idxstats,stats}" into ch_sort_bam_stats

		script:
		"""
		samtools sort $unsorted_bam -T ${prefix} -o ${prefix}.sorted.bam 
		samtools index ${prefix}.sorted.bam
		samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
		samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
		samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
		"""
	}
}

/*
 * Merging replicates & filtering
 */

// Merging
process markDuplicates{
	tag "${prefix}"
	publishDir path: "${params.outdir}/marked_bams", mode: 'copy',
				saveAs: {filename ->
							if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "samtools_stats/$filename"
							else if (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) filename
							else null 
						}
	
	input:
		set val(prefix), file(sorted_bams) from ch_sort_bam

	output:
		set val(prefix), file("${prefix}_marked.{bam,bam.bai}") into ch_marked_bams, ch_marked_preseq
		file "${prefix}_marked.{flagstat,idxstats,stats}" into ch_marked_bam_samstats
		file "*.txt" into ch_marked_bam_picstats

	script:
    bam_files = sorted_bams.findAll { it.toString().endsWith('.bam') }.sort()
	"""
	picard -Xmx4g MarkDuplicates \\
		INPUT=${bam_files[0]} \\
		OUTPUT=${prefix}_marked.bam \\
		ASSUME_SORTED=true \\
		REMOVE_DUPLICATES=false \\
		METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
		VALIDATION_STRINGENCY=LENIENT \\
		TMP_DIR=tmpstararked.bam
	samtools idxstats ${prefix}_marked.bam > ${prefix}_marked.idxstats
	samtools flagstat ${prefix}_marked.bam > ${prefix}_marked.flagstat
	samtools stats ${prefix}_marked.bam > ${prefix}_marked.stats
	"""
}

// Preseq complexity analysis before filtering
if(!params.skip_preseq) {
	process preseqAnalysis{
		tag "${prefix}"
		publishDir "${params.outdir}/preseq"

		input:
		set val(prefix), file(bam) from ch_marked_preseq

		output:
		file "*.ccurve.txt" into ch_preseq_stats

		script:
		"""
		preseq lc_extrap -v -output ${prefix}.ccurve.txt -bam ${bam[0]}
		"""
	}
}

// Filtering
if (!params.skip_filtering){
	process bamFiltering{
		tag "${prefix}"
		publishDir path: "${params.outdir}/filtered_bams", mode: 'copy',
					saveAs: {filename ->
							if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "markstep_stats/$filename"
							else if (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) filename
							else null 
						}
		
		input:
		set val(prefix), file(marked_bam) from ch_marked_bams
		file bed from ch_genome_filter_regions.collect()
		file bamtools_filter_config from ch_bamtools_filter_config

		output:
		set val(prefix), file("*.{bam,bam.bai}") into ch_filtered_bams
		set val(prefix), file("*.flagstat") into ch_filtered_flagstat
		file "*.{idxstats,stats}" into ch_filtered_stats

		script:
		filter_params = params.singleEnd ? "-F 0x004" : "-F 0x004 -F 0x0008 -f 0x001"
		dup_params = params.keep_dups ? "" : "-F 0x0400"
		multimap_params = params.keep_multi_map ? "" : "-q 1"
		blacklist_params = params.blacklist ? "-L $bed" : ""
		name_sort_bam = params.singleEnd ? "" : "samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $prefix ${prefix}_filtered.bam"
		"""
		samtools view \\
			$filter_params \\
			$dup_params \\
			$multimap_params \\
			$blacklist_params \\
			-b ${marked_bam[0]} \\
			| bamtools filter \\
				-out ${prefix}_filtered.bam \\
				-script $bamtools_filter_config
		samtools index ${prefix}_filtered.bam
		samtools flagstat ${prefix}_filtered.bam > ${prefix}_filtered.bam.flagstat
		samtools idxstats ${prefix}_filtered.bam > ${prefix}_filtered.bam.idxstats
		samtools stats ${prefix}_filtered.bam > ${prefix}_filtered.bam.stats
		$name_sort_bam
		"""
	}
}

// Remove orphans from PE files
if (params.singleEnd) {
    ch_filtered_bams
        .into { ch_filtered_bams_metrics;
                ch_filtered_bams_bigwig;
                ch_filtered_bams_macs_1;
                ch_filtered_bams_macs_2;
                ch_filtered_bams_phantompeakqualtools;
                ch_filtered_bams_counts }

    ch_filtered_flagstat
        .into { ch_filtered_flagstat_bigwig;
                ch_filtered_flagstat_macs;
                ch_filtered_flagstat_mqc }

    ch_filtered_stats
        .set { ch_filtered_stats_mqc }
} else {
    process removeOrphanPEbam {
        tag "$prefix"
        publishDir path: "${params.outdir}/rm_orphanPEbam", mode: 'copy',
					saveAs: {filename ->
								if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "samtools_stats/$filename"
								else if (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) filename
								else null 
							}

        input:
        set val(prefix), file(filtered_bams) from ch_filtered_bams

        output:
        set val(prefix), file("*_rm_orphanPEbam.{bam,bam.bai}") into ch_filtered_bams_metrics,
                                                           ch_filtered_bams_bigwig,
                                                           ch_filtered_bams_macs_1,
                                                           ch_filtered_bams_macs_2,
                                                           ch_filtered_bams_phantompeakqualtools
        set val(prefix), file("${prefix}.bam") into ch_filtered_bams_counts
        set val(prefix), file("*.flagstat") into ch_filtered_flagstat_bigwig,
                                               ch_filtered_flagstat_macs,
                                               ch_filtered_flagstat_mqc
        file "*.{idxstats,stats}" into ch_filtered_stats_mqc

        script:
        """
        bampe_rm_orphan.py ${filtered_bams[0]} ${prefix}.bam --only_fr_pairs
        samtools sort -o ${prefix}_rm_orphanPEbam.bam -T $prefix ${prefix}.bam
        samtools index ${prefix}_rm_orphanPEbam.bam
        samtools flagstat ${prefix}_rm_orphanPEbam.bam > ${prefix}_rm_orphanPEbam.bam.flagstat
        samtools idxstats ${prefix}_rm_orphanPEbam.bam > ${prefix}_rm_orphanPEbam.bam.idxstats
        samtools stats ${prefix}_rm_orphanPEbam.bam > ${prefix}_rm_orphanPEbam.bam.stats
        """
    }
}

/*
 * MultiQC
 */

// Retrieve software from environment
process getSoftwareVersionsPython3{
	when:
  		!params.skip_multiqc

	output:
		file 'soft_versions.txt' into software_versions_yaml

	script:
		"""
		echo $workflow.manifest.version &> soft_versions.txt
		echo $workflow.nextflow.version >> soft_versions.txt
		fastqc --version >> soft_versions.txt
		multiqc --version >> soft_versions.txt
		echo \$(bwa 2>&1) >> soft_versions.txt
		bowtie2 --version >> soft_versions_txt
		samtools --version >> soft_versions.txt
		STAR --version >> soft_versions.txt
		echo \$(picard MarkDuplicates --version 2>&1) >> soft_versions.txt
		preseq --version >> soft_versions.txt
		python -c "import pysam; print(pysam.__version__)" >> soft_versions.txt
		scrape_software_versions.py &> software_versions_mqc.yaml
		"""
}

// Workflow summary
process workflow_summary_mqc {
	when:
		!params.skip_multiqc

	output:
	file 'workflow_summary_mqc.yaml' into workflow_summary_yaml

  	exec:
def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
yaml_file.text  = """
id: 'summary'
description: " - this information is collected when the pipeline is started."
section_name: 'Workflow Summary'
section_href: 'https://gitlab.curie.fr/chipseq'
plot_type: 'html'
data: |
	<dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
	</dl>
""".stripIndent()
}


process multiqc {
	publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
		file splan from ch_splan.collect()
		file metadata from ch_metadata.ifEmpty([])
		file multiqc_config from ch_multiqc_config
		file ('software_versions/*') from software_versions_yaml.collect()
		file ('workflow_summary/*') from workflow_summary_yaml.collect()

		file ('fastQC/*') from ch_fastqc_results.collect().ifEmpty([])

		file ('sorted_bams/samtools_stats/*') from ch_sort_bam_stats.collect().ifEmpty([])
		file ('marked_bams/samtools_stats/*') from ch_marked_bam_samstats.collect().ifEmpty([])
		file ('marked_bams/samtools_stats/*') from ch_marked_bam_picstats.collect().ifEmpty([])
		file ('rm_orphanPEbam/samtools_stats/*') from ch_filtered_flagstat_mqc.collect().ifEmpty([])
		file ('rm_orphanPEbam/samtools_stats/*') from ch_filtered_stats_mqc.collect().ifEmpty([])
		
		file ('preseq/*') from ch_preseq_stats.collect().ifEmpty([])

    output:
		file splan
		file "*multiqc_report.html" into multiqc_report
		file "*_data"

    script:
	rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
	rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
	metadata_opts = params.metadata ? "--metadata ${metadata}" : ""
	splan_opts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
	modules_list = "-m custom_content -m fastqc -m samtools -m picard -m preseq"

	"""
	mqc_header.py --name "Chip-seq" --version ${workflow.manifest.version} ${metadata_opts} > multiqc-config-header.yaml
	multiqc . -f $rtitle $rfilename -c $multiqc_config -c multiqc-config-header.yaml $modules_list
	"""
}

/* Creates a file at the end of workflow execution */
workflow.onComplete {
    /*pipeline_report.html*/

    def report_fields = [:]
    report_fields['version'] = workflow.manifest.version
    report_fields['runName'] = custom_runName ?: workflow.runName
    report_fields['success'] = workflow.success
    report_fields['dateComplete'] = workflow.complete
    report_fields['duration'] = workflow.duration
    report_fields['exitStatus'] = workflow.exitStatus
    report_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    report_fields['errorReport'] = (workflow.errorReport ?: 'None')
    report_fields['commandLine'] = workflow.commandLine
    report_fields['projectDir'] = workflow.projectDir
    report_fields['summary'] = summary
    report_fields['summary']['Date Started'] = workflow.start
    report_fields['summary']['Date Completed'] = workflow.complete
    report_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    report_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) report_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) report_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) report_fields['summary']['Pipeline Git branch/tag'] = workflow.revision

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/oncomplete_template.txt")
    def txt_template = engine.createTemplate(tf).make(report_fields)
    def report_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/oncomplete_template.html")
    def html_template = engine.createTemplate(hf).make(report_fields)
    def report_html = html_template.toString()

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << report_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << report_txt }

    /*oncomplete file*/

    File woc = new File("${params.outdir}/workflow.oncomplete.txt")
    Map endSummary = [:]
    endSummary['Completed on'] = workflow.complete
    endSummary['Duration']     = workflow.duration
    endSummary['Success']      = workflow.success
    endSummary['exit status']  = workflow.exitStatus
    endSummary['Error report'] = workflow.errorReport ?: '-'

    String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
    println endWfSummary
    String execInfo = "Execution summary\n${logSep}\n${endWfSummary}\n${logSep}\n"
    woc.write(execInfo)

    /*final logs*/
    if(workflow.success){
        log.info "[Chip-seq] Pipeline Complete"
    }else{
        log.info "[Chip-seq] FAILED: $workflow.runName"
    }
}

import colorlog.MonoChrome
import colorlog.PolyChrome
import groovy.text.GStringTemplateEngine
import nextflow.script.BaseScript
import nextflow.script.WorkflowMetadata
import ParamsLinter
import ParamsReader
import org.slf4j.LoggerFactory
import org.slf4j.Logger
import nextflow.Nextflow
import nextflow.Channel


class NFTools {

    static final Logger log = LoggerFactory.getLogger(Nextflow.class)

    private static LinkedHashMap generateLogColors(Boolean monochromeLogs = false) {
         Map colors = monochromeLogs ? new MonoChrome().palette() : new PolyChrome().palette()
         return colors as LinkedHashMap
     }


    /**********************************************
     *
     * Welcome and Logs
     *
     **********************************************/

    /**
     * Welcome method which should be launch at the beginning of the workflow
     * @return
     */

    public static void welcome(workflow, params) {
        def colors = generateLogColors(params.get("monochromeLogs", false) as Boolean)
        nfHeader(params, workflow as WorkflowMetadata)
        if ("${workflow.manifest.version}" =~ /dev/ ) {
	   //def devMessageFile = new File("${workflow.projectDir}/assets/devMessage.txt")
           //log.info devMessageFile.text
           printDisclaimer(params, workflow as WorkflowMetadata)
        }
        if (params.help) {
            def paramsWithUsage = readParamsFromJsonSettings("${workflow.projectDir}/parameters.settings.json")
            helpMessage(paramsWithUsage, workflow)
            Nextflow.exit(1)
        }
    }

    /**
     * Print Header in assets in the log
     * @param params
     * @param workflow
     * @return
     */

    public static void nfHeader(params, WorkflowMetadata workflow) {

        LinkedHashMap context = generateLogColors(params.get("monochromeLogs", false) as Boolean)
        context << workflow.properties

        def engine = new GStringTemplateEngine()
        def txtTemplate = engine.createTemplate(
          new File("${workflow.projectDir}/assets/nfHeader.txt")
        ).make(context)
	log.info txtTemplate.toString()
    }

    /**
     * Print Header in assets in the log
     * @param params
     * @param workflow
     * @return
     */

    public static void printDisclaimer(params, WorkflowMetadata workflow) {

        LinkedHashMap context = generateLogColors(params.get("monochromeLogs", false) as Boolean)
        context << workflow.properties

        def engine = new GStringTemplateEngine()
        def txtTemplate = engine.createTemplate(
          new File("${workflow.projectDir}/assets/devMessage.txt")
        ).make(context)
        log.info txtTemplate.toString()
    }



    /**
     * Has the run name been specified by the user?
     * This has the bonus effect of catching both -name and --name
     *
     * @return customRunName
     */
    
     public static String checkRunName(workflowRunName, runName) {
       return workflowRunName ==~ /[a-z]+_[a-z]+/ && runName ? runName : workflowRunName
     }


    /*************************************************
     *
     * Help functions
     *
     *************************************************/

     private static Map formatParameterHelpData(params) {

        Map result = [name: params.get('name'), value: '', usage: params.get('usage')]
        // value describes the expected input for the param
        result.value =  [params.type == boolean.toString() ? '' : params.type.toString().toUpperCase(), params.choices ?: ''].join(' ')
        return result
     }

     private static String prettyFormatParamGroupWithPaddingAndIndent(List paramGroup, String groupName, Integer padding = 2, Integer indent = 4) {

        def maxParamNameLength = paramGroup.collect { it.name.size() }.max()
        def paramChoices = paramGroup.findAll { it.choices }.collect { it.choices }
        def maxChoiceStringLength = paramChoices.collect { it.toString().size() }.max() ?: 0
        def maxTypeLength = paramGroup.collect { (it.type as String).size() }.max() ?: 0
        def maxValueLength = maxChoiceStringLength + maxTypeLength
        def usagePadding = indent + maxParamNameLength + maxValueLength + 2 * padding + 2
        Integer usageLength = ((150 - usagePadding) / 100).round(1) * 100

        def paramsFormattedList = paramGroup.sort { it.name }.collect {
            Map param ->
                def paramHelpData = formatParameterHelpData(param)
                def usage = []

                "${paramHelpData.usage}".eachMatch(".{1,${usageLength - 1}}([\\s\\.]|\$)", {
                    usage << it.pop()
                })
                usage = usage.join("\n" + " " * usagePadding)
                sprintf("%${indent}s%-${maxParamNameLength + padding}s %-${maxValueLength + padding}s %s\n", "", "--${paramHelpData.name}", "${paramHelpData.value}", "${usage}")

        }
        String.format("%s:\n%s", groupName.toUpperCase(), paramsFormattedList.join()).stripIndent()
     }

     private static String prettyFormatParamsWithPaddingAndIndent(List paramsWithUsage, Integer padding = 2, Integer indent = 4) {

        def groupedParamsWithUsage = paramsWithUsage.groupBy { it.group }
        def formattedParamsGroups = groupedParamsWithUsage.collect {
            prettyFormatParamGroupWithPaddingAndIndent(it.value, it.key, padding, indent)
        }
        return formattedParamsGroups.join('\n')
     }


    /**
     * Print Profile Help message
     */


    static String profileHelp(){

    String.format("""\
=======================================================
Available Profiles
   -profile test                        Run the test dataset
   -profile conda                       Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
   -profile multiconda                  Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
   -profile path                        Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
   -profile multipath                   Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path
   -profile docker                      Use the Docker images for each process
   -profile singularity                 Use the Singularity images for each process. Use `--singularityPath` to define the insallation path
   -profile cluster                     Run the workflow on the cluster, instead of locally
""")
    }


    /**
     * Generate Help message
     * @param paramsWithUsage
     * @param workflow
     * @return
     */

    static String helpMessage(paramsWithUsage, workflow) {
        def CLIHelpMsg = []
        paramsWithUsage.each {it.group == "Mandatory arguments" ? CLIHelpMsg << "--" + it.name << it.type.toUpperCase() : ""}
        log.info String.format("""\

            Usage:

            The typical command for running the pipeline is as follows:

            nextflow run main.nf ${CLIHelpMsg.join(" ")}

        %s
        %s
            """.stripIndent(), prettyFormatParamsWithPaddingAndIndent(paramsWithUsage, 2, 4), profileHelp())
    }



    /***********************************
     *
     * Checks params with json
     *
     **********************************/

    /**
     * Lint params scope with usage from the nf-core json
     * @param params
     * @param paramsWithUsage
     * @return
     */

    public static Object lint(params, paramsWithUsage) {
    	      def linter = new ParamsLinter(params, paramsWithUsage, log)
    	         return linter.lint()
    }

    /**
     * Reads params in json file
     * @param path
     * @return
     */

    static Object readParamsFromJsonSettings(String path) {

        def paramsWithUsage
        try {
            paramsWithUsage = ParamsReader.readParamsFromJsonSettings(path).get("parameters")
        } catch (Exception e) {
            log.warn "Could not read parameters settings from JSON. $e"
            paramsWithUsage = Collections.emptyMap()
	        }
        return paramsWithUsage
    }


    /****************************************
     *
     * Load configuration from genome.config
     *
     ***************************************/

    /**
     * Get variables from genome config file
     * @param params
     * @param value
     * @return
     */

     public static String getGenomeAttribute(params, attribute, genome=params.genome) {
        def val = ''
        if (params.genomes && params.genomes.containsKey(genome)) {
            if (params.genomes[ genome ].containsKey(attribute)) {
                val = params.genomes[ genome ][ attribute ]
            }
        }
        return val
    }



    /*****************************************
     *
     * Summary
     *
     ****************************************/
    /**
     * Print summary dict with the logger and return formatted channel
     * @param summary
     * @param workflow
     * @return
     */

     public static Object summarize(summary, workflow, params) {
        def colors = generateLogColors(params.get("monochromeLogs", false) as Boolean)
        log.info summary.collect { k, v -> "${k.padRight(18)}: $v" }.join("\n") +
                "\n${colors.dim}" + "-"*80 + "${colors.reset}"
        return Channel.fromList(summary.collect{ [it.key, it.value] })
          .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
                .reduce { a, b -> return [a, b].join("\n            ") }
                .map { x -> """
    id: 'summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'Workflow Summary'
    section_href: '$workflow.manifest.homePage'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent()}
     }



     /****************************************
      *
      * Load input data
      *
      ****************************************/

     /**
      * Check if path has the correct file extension
      * @param it
      * @param extension
      * @return
      */

      public static boolean hasExtension(it, String extension) {
        it.toString().toLowerCase().endsWith(extension.toLowerCase())
      }

      /**
       * Check if file extension from input path is defined in allowed extensions list
       * @param path
       * @param extensions
       * @return
       */

      public static String getExtension(path, List extensions, params) {
        for (extension in extensions) {
            if (hasExtension(path, extension)) {
                return extension
            }
        }

        def colors = generateLogColors(params.get("monochromeLogs", false) as Boolean)
        log.warn """\
        Can't guess field separator with the actual input file extension. Trying with the default one (.csv)
        """.stripIndent().toString()
        return ""
      }

      /**
       * Check if a row has the expected number of item
       *
       * @param row
       * @param number
       * @return Boolean
       */

      public static boolean checkNumberOfItem(row, Integer number, params) {
        def colors = generateLogColors(params.get("monochromeLogs", false) as Boolean)
        if (row.size() != number) {
           def message = """\
               ${colors.red}[WARNING] Malformed row in input file: ${row}
               Input file should have ${number} but have ${row.size()} items. see --help for more information${colors.reset}
               """.stripIndent()
           log.info message.toString()
        }
        return true
      }

      /**
       * Return file if it exists
       *
       * @return Nextflow.file Object
       */
      public static Object returnFile(String it, params) {
        def colors = generateLogColors(params.get("monochromeLogs", false) as Boolean)
        if (it =~ /(http|ftp)/) {
           return Nextflow.file(it)
        } else if (!Nextflow.file(it).exists()) {
           Nextflow.exit(
             "${colors.red}[WARNING] Input file does not exists: ${it}, see --help for more information${colors.reset}"
           )
        } else {
           return Nextflow.file(it)
        }
      }

      /**
       * Channeling the input file containing FASTQ or BAM
       * Format is: "idSample,sampleName,pathToFastq1,[pathToFastq2]"
       * or: "sampleID,sampleName,pathToBam"
       *
       * @param samplePlan
       * @return
       */

      public static Object getInputData(samplePlan, reads, readPaths, singleEnd, params) {

        if (samplePlan) {
      	  return Channel
            .fromPath(samplePlan)
            .splitCsv(header: false)
            .map { row ->
              def sampleID = row[0]
              def sampleName = row[1]
              def inputFile1 = returnFile(row[2], params)
              def inputFile2 = 'null'

              if ((!singleEnd) && (hasExtension(inputFile1, 'fastq.gz') || hasExtension(inputFile1, 'fq.gz') || hasExtension(inputFile1, 'fastq'))) {
                checkNumberOfItem(row, 4, params)
                inputFile2 = returnFile(row[3], params)
                if (!hasExtension(inputFile2, 'fastq.gz') && !hasExtension(inputFile2, 'fq.gz') && !hasExtension(inputFile2, 'fastq')) {
                  Nextflow.exit(1, "File: ${inputFile2} has the wrong extension. See --help for more information")
                }
              } else if (hasExtension(inputFile1, 'bam')) {
                checkNumberOfItem(row, 3, params)
              } else {
                log.warn "No recognisable extention for input file: ${inputFile1}"
              }
              return singleEnd ? [sampleID, [inputFile1]] : [sampleID, [inputFile1, inputFile2]]
            }
        } else if (readPaths) {
          return Channel
            .fromList(readPaths)
            .map { row ->
              def sampleId = row[0]
              def inputFile1 = returnFile(row[1][0], params)
              def inputFile2 = singleEnd ? null: returnFile(row[1][1], params)
              singleEnd ? [sampleId, [inputFile1]] : [sampleId, [inputFile1, inputFile2]]
           }.ifEmpty { Nextflow.exit 1, "params.readPaths was empty - no input files supplied" }
        } else {
          return Channel
            .fromFilePairs(reads, size: singleEnd ? 1 : 2)
            .ifEmpty { Nextflow.exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
            .map { row -> singleEnd ? [row[0], [row[1][0]]] : [row[0], [row[1][0], row[1][1]]] }
        }
      }


      /**
       * Channeling the samplePlan and create a file is no samplePlan is provided
       * 
       * @param samplePlan
       * @param reads
       * @param readPaths
       * @param singleEnd
       
       * @return
       */

      public static Object getSamplePlan(samplePlan, reads, readPaths, singleEnd) {
        if (samplePlan){
	  return Channel.fromPath(samplePlan)
	} else if(readPaths){
          if (singleEnd){
	    return Channel
	      .from(readPaths)
              .collectFile() {
	        item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
               }
          }else{
            return Channel
              .from(readPaths)
              .collectFile() {
                item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
              }
          }
        }else{
	  if (singleEnd){
	    return Channel
	      .fromFilePairs( reads, size: 1 )
	      .collectFile() {
	        item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
	      }
  	  }else{
	    return Channel
	      .fromFilePairs( reads, size: 2 )
	      .collectFile() {
	        item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
 	      }
	  }
        }
      }


   /************************************
    *
    * Workflow on complete
    *
    ************************************/

    /**
     * Generate reports at the end of the workflow
     * @param workflow
     * @param params
     * @param reportFields
     * @param multiQCOutCh
     * @return
     */

    public static void makeReports(workflow, params, summary, customRunName, mqcReport) {
    	def colors = generateLogColors(params.get("monochromeLogs", false) as Boolean)  
        def engine = new groovy.text.GStringTemplateEngine()

	// Complete summary
	summary: summary + [
      	'Date Started': workflow.start,
	'Date Completed': workflow.complete,
	'Pipeline script file path': workflow.scriptFile,
	'Pipeline script hash ID': workflow.scriptId,
	'Pipeline repository Git URL': workflow.repository ?: null,
	'Pipeline repository Git Commit': workflow.commitId ?: null,
	'Pipeline Git branch/tag': workflow.revision ?: null,
	'Docker image': workflow.container ?: null,
	'Nextflow Version': workflow.nextflow.version,
	'Nextflow Build': workflow.nextflow.build,
	'Nextflow Compile Timestamp': workflow.nextflow.timestamp,
	].findAll{ it.value != null }

  	def reportFields = params + [
        workflowName: workflow.manifest.name,
	version : workflow.manifest.version,
        runName: customRunName, 
 	success: workflow.success,
	dateComplete: workflow.complete,
	duration:workflow.duration,
	exitStatus: workflow.exitStatus,
	errorMessage: (workflow.errorMessage ?: 'None'),
	errorReport: (workflow.errorReport ?: 'None'),
	commandLine: workflow.commandLine,
	projectDir: workflow.projectDir,
	summary: summary]

	// Render the TXT template
        def txtTemplate = engine.createTemplate(new File("$workflow.projectDir/assets/onCompleteTemplate.txt")).make(reportFields)
        def txtReport = txtTemplate.toString()

        // Render the HTML template
        def hf = new File("$workflow.projectDir/assets/onCompleteTemplate.html")
        def htmlTemplate = engine.createTemplate(hf).make(reportFields)
        def htmlReport = htmlTemplate.toString()

        // Write summary e-mail HTML to a file
        def summaryDir = new File("${params.summaryDir}")
        if (!summaryDir.exists()) summaryDir.mkdirs()
        def output_hf = new File(summaryDir, "pipelineReport.html")
        output_hf.withWriter { w -> w << htmlReport }
        def output_tf = new File(summaryDir, "pipelineReport.txt")
        output_tf.withWriter { w -> w << txtReport }

        // On success try attach the multiqc report
        try {
            if (workflow.success) {
                if (mqcReport.getClass() == ArrayList) {
                    //log.warn "[$workflow.manifest.name] Found multiple reports from process 'multiqc', will use only one"
                    mqcReport = mqcReport[0]
                }
            }
        } catch (all) {
            log.warn "[$workflow.manifest.name] Could not attach MultiQC report to summary email"
        }

        if (params.email) {
            // Set up the e-mail variables
            def subject = workflow.success? "[$workflow.manifest.name] Successful: ${workflow.runName}": "[$workflow.manifest.name] FAILED: ${workflow.runName}"
            // Render the sendmail template
            def smailFields = [
                    to: params.email,
                    subject: subject,
                    text: txtReport,
                    body: htmlReport,
                    attach: mqcReport,
            ]
            Nextflow.sendMail(smailFields)
            log.info "[$workflow.manifest.name] Sent summary e-mail to $params.email (sendmail)"
        }

        // workflowOnComplete file
 	def outDir = new File("${params.outDir}")
        if (!outDir.exists()) outDir.mkdirs()
        File woc = new File(outDir, "workflowOnComplete.txt")
        def endSummary = [
            'Completed on': workflow.complete,
            'Duration': workflow.duration,
            'Success': workflow.success,
            'Exit status': workflow.exitStatus,
            'Error report': workflow.errorReport ?: '-'
        ]
        
        String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
        String execInfo = "Execution summary\n${endWfSummary}\n"
        woc.withWriter { w -> w << execInfo }

        def endMessage = workflow.success ? workflow.stats.ignoredCount > 0 ? """\
            ${colors.purple}Warning, pipeline completed, but with errored process(es)${colors.reset}
            ${colors.red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt}${colors.reset}
            ${colors.green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt}${colors.reset}
            [$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}
            """.stripIndent() : "[$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}": "[$workflow.manifest.name]${colors.red} Pipeline completed with errors${colors.reset}"
        log.info endMessage.toString()
    }

}
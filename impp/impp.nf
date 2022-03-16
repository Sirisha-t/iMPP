
  #!/usr/bin/env nextflow

/**
integrated Metagenomic Protein Predictor (iMPP)

*/


def versionMessage()
{
	log.info"""

	iMPP - Version: ${workflow.manifest.version}
	""".stripIndent()
}

def helpMessage()
{
	log.info"""

  Usage:
  nextflow run impp.nf --1 R1.fastq --2 R2.fastq --prefix prefix --outdir path [options]

  Input data options:
   --12            <filename>    : fastq file with interlaced forward and reverse paired-end reads
   --1              <filename>    : fastq file with forward paired-end reads
   --2             <filename>    : fastq file with reverse paired-end reads
   --s             <filename>    : fastq file with unpaired reads
   --prefix        <name>         : prefix used to name the results file
   --o/--outdir     <dirname>     : [required] output directory name including the full path
   --n              <int>         : [optional] number of threads
   --m/--max-len    <int>         : [optional] maximum extension length for anchors [default: 150]
   --p/--param-file <filename>    : [optional] parameter file
   --h/--help                     :  help message

  NOTE: Input FASTQ read file must be specified in either --12, --1 & --2 or --s format.


iMPP supports FASTQ and compressed FASTQ files.
"""
}

/**
Prints version when asked for
*/
if (params.version) {
	versionMessage()
	exit 0
}

/**
Prints help when asked for
*/

if (params.help) {
	helpMessage()
	exit 0
}

/**
STEP 0.

Checks input parameters and (if it does not exists) creates the directory
where the results will be stored (aka working directory).
Initialises the log file.

The working directory is named after the prefix and located in the outdir
folder. The log file, that will save summary statistics, execution time,
and warnings generated during the pipeline execution, will be saved in the
working directory as "prefix.log".
*/


//Checking user-defined parameters
if ( !params.single_file && !params.forward_file && !params.$reverse_file && !params.$interleaved_file) {
	exit 1, "Input parameters were not specified properly. Please specify FASTQ reads as input in either single (--s) or paired mode (--1, --2 or --12)."
}


//--reads2 can be omitted when the library layout is "single" (indeed it specifies single-end
//sequencing)
if ( !params.single && (params.reads2 == "null") ) {
	exit 1, "If dealing with paired-end reads, please set the reads2 parameters\nif dealing with single-end reads, please set the library layout to 'single'"
}

//--reads1 and --reads2 can be omitted (and the default from the config file used instead)
//only when mode is "characterisation". Obviously, --reads2 should be always omitted when the
//library layout is single.
if (( (!params.singleEnd && (params.reads1 == "null" || params.reads2 == "null")) || (params.singleEnd && params.reads1 == "null")) ) {
	exit 1, "Please set the reads1 and/or reads2 parameters"
}

//Creates working dir
workingpath = params.outdir + "/" + params.prefix
workingdir = file(workingpath)
if( !workingdir.exists() ) {
	if( !workingdir.mkdirs() ) 	{
		exit 1, "Cannot create working directory: $workingpath"
	}
}


// Header log info
log.info """---------------------------------------------
---------------------------------------------
Analysis introspection:
"""

def summary = [:]

summary['Starting time'] = new java.util.Date()
//Environment
summary['Environment'] = ""
summary['Pipeline Name'] = 'iMPP'
summary['Pipeline Version'] = workflow.manifest.version

summary['Config Profile'] = workflow.profile
summary['Resumed'] = workflow.resume

summary['Nextflow version'] = nextflow.version.toString() + " build " + nextflow.build.toString() + " (" + nextflow.timestamp + ")"

summary['Java version'] = System.getProperty("java.version")
summary['Java Virtual Machine'] = System.getProperty("java.vm.name") + "(" + System.getProperty("java.vm.version") + ")"

summary['Operating system'] = System.getProperty("os.name") + " " + System.getProperty("os.arch") + " v" +  System.getProperty("os.version")
summary['User name'] = System.getProperty("user.name") //User's account name

summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container


//General
summary['Running parameters'] = ""
summary['Reads'] = "[" + params.reads1 + ", " + params.reads2 + "]"
summary['Prefix'] = params.prefix
summary['Running mode'] = params.mode
summary['Layout'] = params.singleEnd ? 'Single-End' : 'Paired-End'


log.info summary.collect { k,v -> "${k.padRight(27)}: $v" }.join("\n")
log.info ""


/**
	Prepare workflow introspection
	This process adds the workflow introspection (also printed at runtime) in the logs
	This is NF-CORE code.
*/

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'workflow-summary'
    description: "This information is collected when the pipeline is started."
    section_name: 'IMPP Workflow Summary'
    section_href: 'https://github.com/alesssia/yamp'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd>$v</dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


if (params.singleEnd) {
	Channel
	.from([[params.prefix, [file(params.reads1)]]])
	.into { read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants }
} else {
	Channel
	.from([[params.prefix, [file(params.reads1), file(params.reads2)]]] )
	.into { read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants }
}

process runimpp{

	input:
	tuple val(name), file(reads) from read_files_dedup

	output:
	tuple val(name), path("${name}_dedup*.fq.gz") into to_synthetic_contaminants

	when:
	params.mode != "characterisation" && params.dedup

	script:
	// This is to deal with single and paired end reads
	def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
	def output = params.singleEnd ? "out=\"${name}_dedup.fq.gz\"" :  "out1=\"${name}_dedup_R1.fq.gz\" out2=\"${name}_dedup_R2.fq.gz\""

	"""
    perl impp_run.pl -Xmx\"\$maxmem\" --s $input --o $output

	"""
}


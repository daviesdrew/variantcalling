#!/usr/bin/env nextflow

nextflow.preview.dsl=2
/*
========================================================================================
                         nf-core/illuminavariantcalling
========================================================================================
 nf-core/illuminavariantcalling Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/illuminavariantcalling
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    
    //colour scheme parameters
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_bold = params.monochrome_logs ? '' : "\033[1m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_block = params.monochrome_logs ? '' : "\033[3m";
    c_ul = params.monochrome_logs ? '' : "\033[4m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_bul = c_bold + c_ul;

    //log.info nfcoreHeader()
    log.info"""
    =${c_dim}===============================================================${c_reset}
    ${c_blue+c_bold}${workflow.manifest.name}${c_reset}   ~  version ${c_purple}${workflow.manifest.version}${c_reset} 
    ${c_dim}================================================================${c_reset}

    ${c_ul}Git info:${c_reset} $workflow.repository - $workflow.revision [$workflow.commitId]

    ${c_bul}Usage:${c_reset}
    The typical command for running the pipeline is as follows:
        
        nextflow run ${workflow.manifest.name}

    The above assumes that you have minimap2 and bowtie2 installed and can
    be located in the path of your environment

    It further assumes that your environment has Freebayes, bcftools and 
    SnpEff installed. 

    nextflow run nf-core/illuminavariantcalling --reads '*_R{1,2}.fastq.gz' -profile docker

    ${c_bul}Mandatory Options:${c_reset}
      --reads [file]                Path to input data (must be surrounded with quotes)
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

//=============================================================================
// HELPER FUNCTIONS
//=============================================================================

def check_dict_args(tool, args, acceptable_args) {
    args.keySet()
        .each{ arg -> assert acceptable_args.contains(arg) }
    
    return args
}

def arg_str_to_dict(arg_str, dict) {
    arg_str.tokenize(',')
        .each{ arg -> 
            dict[arg.tokenize()[0].replace('-','')] = arg.tokenize()[1] }
}

def print_arg_comparison(tool, args, acceptable_args) {
    println("${tool} selected\n")
    println("${tool} args")
    println("${args}\n")
    println("${tool} accepted args")
    println("${acceptable_args}\n")
}

//*****************************************************************************
// HELPER FUNCTIONS: ALIGNMENT
//*****************************************************************************

//check for installation of aligner
//check for argument validitiy
def check_aligner(aligner) {
    println("Checking aligner arguments")

    align_args = [:]
    arg_str_to_dict(params.align_args, align_args)

    bwa_accept_args = ['key', 'bwa', 'one']
    minimap2_accept_args = ['key', 'minimap2']
    bowtie2_accept_args = ['key', 'bowtie2']
    
    switch(aligner) {
        case "bwa": 
            print_arg_comparison(aligner, 
                                 align_args, 
                                 bwa_accept_args)
            return check_dict_args(aligner, 
                                      align_args, 
                                      bwa_accept_args)
            break;

        case "minimap2":
            print_arg_comparison(aligner, 
                                 align_args, 
                                 minimap2_accept_args)
            return check_dict_args(aligner, 
                                      align_args, 
                                      minimap2_accept_args)
            break;

        case "bowtie2": 
            print_arg_comparison(aligner, 
                                 align_args, 
                                 bowtie2_accept_args)
            return check_dict_args(aligner, 
                                      align_args, 
                                      bowtie2_accept_args)
            break;

        default: 
            print_arg_comparison(aligner, 
                                 align_args, 
                                 bowtie2_accept_args)
            return check_dict_args(aligner, 
                                      align_args, 
                                      bowtie2_accept_args)
    }
}

//*****************************************************************************

//*****************************************************************************
// HELPER FUNCTIONS: VARIANT CALLING
//*****************************************************************************

def check_variant_caller(variant_caller) {
    println("Checking variant caller arguments")
    variant_caller_args = [:]
    arg_str_to_dict(params.variant_args, variant_caller_args)
    
    freebayes_accept_args = ['key', 'freebayes']
    
    print_arg_comparison(variant_caller, 
                         variant_caller_args, 
                         freebayes_accept_args)
    check_dict_args(variant_caller, 
                    variant_caller_args, 
                    freebayes_accept_args)
}



//*****************************************************************************

//=============================================================================
// INPUT VALIDATION
//=============================================================================

ref_aligner = (params.containsKey('align') && 
               params.align != null) ? 
                check_aligner(params.align) : null

var_caller = (params.containsKey('variant') && 
                  params.variant != null) ? 
                  check_variant_caller(params.variant) : null


//=============================================================================
// WORKFLOW RUN PARAMETERS LOGGING
//=============================================================================

custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ) {
    custom_runName = workflow.runName
}

def summary = [:]
summary['Pipeline Name']        = workflow.manifest.name
summary['Pipeline Version']     = workflow.manifest.version
summary['Run Name']             = custom_runName ?: workflow.runName
summary['Reads']                = params.reads

summary['Aligner']              = params.align
summary['Aligner Args']         = params.align_args
summary['Variant Caller']       = params.variant
summary['Variant Caller Args']  = params.variant_args

summary['Max Memory']           = params.max_memory
summary['Max CPUs']             = params.max_cpus
summary['Max Time']             = params.max_time

summary['Container Engine']     = workflow.containerEngine
summary['Container']            = workflow.container

summary['Current Home']         = "$HOME"
summary['Current User']         = "$USER"
summary['Current Path']         = "$PATH"
summary['Working Dir']          = workflow.workDir
summary['Output Dir']           = file(params.outdir)
summary['Script Dir']           = workflow.projectDir
summary['Config Profile']       = workflow.profile

//=============================================================================
// PROCESSES
//=============================================================================

//*****************************************************************************
// PROCESSES: REMOVE_PHIX
// Remove any Coliphage phi-X714 reads using bbduk
// bbduk and options explained above script call
//*****************************************************************************
process REMOVE_PHIX {
    //tag process with sample id
    tag "$sample_id"
    
    //publishes output files making a regexp to 
    //directories described using regexp
    publishDir "${params.outdir}/qc/phix_removal", 
                pattern: "*.txt", mode: 'copy'
    publishDir "${params.outdir}/reads/phix_removed",
                pattern: "*.fastq.qz"

    input: 
        file phix
        tuple sample_id, path(r1), path(r2)

    //the emit argument allows us to access the output of the process
    //using names in dictionary style in the workflow 
    output: 
        tuple sample_id, path(r1_out), path(r2_out), 
                emit: 'reads'
        tuple sample_id, path(stats), 
                emit: 'stats'
   
    //values used as script args
    script:
    r1_out = "${sample_id}_1.phix_removed.fastq.gz"
    r2_out = "${sample_id}_2.phix_removed.fastq.qz"
    stats = "${sample_id}-remove_phix-stats.txt"

    /*
    bbduk (Decontamination using kmers)

    -Xmx${numerical value}(unit size) 
        Defines the maximum size of the heap allocated by the JVM

    -in1 and in2 are input reads 

    -out1 and out2 are output of decontaminated reads

    -ref is the reference genome that the reads are decontaminated of

    -k is kmer length

    -hdist is the maximum allowable hamming distance

    -stats is the output file for statistic analysis
    */ 
    
    """
    bbduk \\
        -Xmx${task.memory.toMega()}m \\
        in1=$r1 in2=$r2 \\
        out1=$r1_out out2=$r2_out \\
        ref=$phix k=31 hdist=1 \\
        stats=$stats
    """
}
//*****************************************************************************

//*****************************************************************************
// PROCESSES: FASTP
// Adapter trimming and quality filering using fastp
// By default, only up to 40% of bases can be below the min Phred score base
// quality threshold of Q15
// fastp and options explained above script call
//*****************************************************************************
process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp/html",
                pattern: "*.html", mode: 'copy'
    publishDir "${params.outdir}/fastp/json",
                pattern: "*.json", mode: 'copy'
    publishDir "${params.outdir}/reads/fastp",
                pattern: "*.fastp.fastq.gz"
                
    input:
        tuple sample_id, path(r1), path(r2)
        
    output: 
        tuple sample_id, path(r1_out), path(r2_out), 
                emit: 'reads'
        tuple sample_id, path(html_report), path(json_report), 
                emit: 'report'
                
    script:
        r1_out = "${sample_id}_1.fastp.fastq.gz"
        r2_out = "${sample_id}_2.fastp.fastq.gz"
        json_report = "fastp-report-${sample_id}.json"
        html_report = "fastp-report-${sample_id}.html"
        
    """
    fastp \\
        -i $r1 -I $r2 \\
        -o $r1_out -O $r2_out \\
        -p -c -R "$sample_id fastp report" \\
        -w ${task.cpus} \\
        -q ${params.fastp_min_base_quality} \\
        -u ${params.fastp_max_percent_low_qual_base} \\
        -j $json_report -h $html_report 
   """
}
//*****************************************************************************

//*****************************************************************************
// PROCESSES: FASTQC
// Quality control using fastqc
// fastqc and options explained above script call
//*****************************************************************************
process FASTQC { 
    tag "$sample_id"
    publishDir "${params.outdir}/qc/fastqc",
                mode: 'copy',
                saveAs: { filename -> 
                    filename.indexOf(".zip") > 0 
                        ? "zips/$filename"
                        : "$filename"
                }

    input: 
        tuple val(sample_id), path(r1), path(r2)

    output:
        file "*_fastqc.{zip,html}"
    
    script:
    """
    fastqc -q $r1 $r2
    """
}
//*****************************************************************************

//*****************************************************************************
// PROCESSES: BWA 
// Paired End read alignment 
// Utilizes reference genome
//*****************************************************************************

process BWA {
    tag "$sample_id"
    publishDir "${params.outdir}/align/bwa", 
                pattern: "*.sam", mode: 'copy'                

    input:
        file ref            
        tuple val(sample_id), path(r1), path(r2)
    
    output:
        tuple val(sample_id), path(r1_out), path(r2_out),
                emit: 'align'
    
    script:
        r1_unzip = "${sample_id}_1.fastp.fastq"
        r2_unzip = "${sample_id}_2.fastp.fastq"
        r1_out = "${sample_id}_1_align.sam"
        r2_out = "${sample_id}_2_align.sam"
        align = "${sample_id}_aln-pe.sam"

    """
    gzip -d --force $r1 > r1_unzip;
    gzip -d --force $r2 > r2_unzip;
    bwa index $ref;
    bwa aln $ref $r1_unzip > $r1_out; 
    bwa aln $ref $r2_unzip > $r2_out;
    bwa sampe $ref $r2_out $r2_out $r1_unzip $r2_unzip > $align
    """
}

//*****************************************************************************

//*****************************************************************************
// PROCESSES: SAMTOBAM
// Convert files from SAM to BAM format  
// Prepares output from BWA to become input for Freebayes
//*****************************************************************************

process SAMTOBAM {
    tag "$sample_id"
    publishDir "${params.outdir}/align/bwa",
                pattern: "*.bam", mode: "copy"

    input: 
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path(r1_out), path(r2_out)

    script: 
        r1_out = "${sample_id}_1.bam"
        r2_out = "${sample_id}_2.bam"

    """
    samtools view -S -b $r1
    """
}

//*****************************************************************************

//=============================================================================
// WORKFLOW DEFINITION 
//=============================================================================

workflow {
    // Create channel for paired end reads
    
    //obtain read files from specified location 
    ch_reads = Channel.fromFilePairs(
        params.reads,
        flat: true,
        checkIfExists: true)
    //check if Channel is empty
    .ifEmpty{ exit 1, "No reads specified! Please specify where your reads are, e.g. '--reads \"/path/to/reads/*R{1,2}.fastq.qz\"' (quotes around reads ath required if using `*` and other characters expanded by the shell!)"}
    //Convert specified files into tuples
    .map{ [ it[0].replaceAll(/_S\d{1,2}_L001/,""), it[1], it[2] ] }
      
    ch_phix = Channel.value(file("${baseDir}/data/phix.fa"))
    
    REMOVE_PHIX(ch_phix, ch_reads)
    FASTP(REMOVE_PHIX.out.reads)
    fastqc_reports = FASTQC(FASTP.out.reads)

    align_ref = Channel.value(file("${baseDir}/data/ref.fa"))

    BWA(align_ref, FASTP.out.reads)
    println("${BWA.out.align}") 
}

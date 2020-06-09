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

//----------------------------------------
// HELPER FUNCTIONS: GENERAL TOOLS
//----------------------------------------
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

def check_args(tool, tool_args, tool_accept_args) {
    print_arg_comparison(tool, tool_args, tool_accept_args)
    return check_dict_args(tool, tool_args, tool_accept_args)
}


//----------------------------------------
// HELPER FUNCTIONS: INPUT VALIDATION
//----------------------------------------

def input_validation(tool_type, tool, tool_args) {
    if (params[tool_type] == null || params[tool_args] == null ||
        !(params.containsKey(tool_type) && params.containsKey(tool_args))) {
        params[tool_type] = null
        params[tool_args] = null
    } else {
        switch(tool_type) {
            
            case 'align':
                check_aligner(tool)
                break;

            case 'variant': 
                check_variant_caller(tool)
                break;

            case 'filter':
                check_filter(tool)
                break;
        
            case 'consensus':
                check_consensus(tool)
                break;
        }
    }
}
//----------------------------------------

//----------------------------------------
// HELPER FUNCTIONS: CHECK ALIGNER ARGS
//----------------------------------------
def check_aligner(aligner) {
    println("Checking aligner arguments")

    align_args = [:]
    arg_str_to_dict(params.align_args, align_args)

    bwa_accept_args = ['key', 'bwa', 'one']
    minimap2_accept_args = ['key', 'minimap2']
    bowtie2_accept_args = ['key', 'bowtie2']
    
    switch(aligner) {
        case "bwa": 
            check_args(aligner, 
                       align_args, 
                       bwa_accept_args)
            break;

        case "minimap2":
            check_args(aligner,
                       align_args, 
                       minimap2_accept_args)
            break;

        case "bowtie2": 
            check_args(aligner, 
                       align_args,
                       bowtie2_accept_args)
            break;

        default: 
            check_args(aligner, 
                       align_args,
                       bowtie2_accept_args)
    }
}
//----------------------------------------

//----------------------------------------
// HELPER FUNCTIONS: CHECK VARIANT CALLING ARGS
//----------------------------------------
def check_variant_caller(variant_caller) {
    println("Checking variant caller arguments")
    variant_caller_args = [:]
    arg_str_to_dict(params.variant_args, variant_caller_args)
    
    freebayes_accept_args = ['key', 'freebayes']
    
    switch(variant_caller) {
        case "freebayes": 
            check_args(variant_caller, 
                      variant_caller_args,
                      freebayes_accept_args)
            break;

        default: 
           check_args(variant_caller, 
                      variant_caller_args,
                      freebayes_accept_args)
    } 
}
//----------------------------------------

//----------------------------------------
// HELPER FUNCTIONS: CHECK CONSENSUS ARGS
//----------------------------------------
def check_consensus(consensus) {
    println("Checking consensus arguments")
    consensus_args = [:]
    arg_str_to_dict(params.consensus_args, consensus_args)
    
    bcftools_accept_args = ['key', 'bcftools']
    vcf_consensus_accept_args = ['key', 'vcf_consensus']

    switch(consensus) {
        case "bcftools": 
            check_args(consensus, 
                       consensus_args,
                       bcftools_accept_args)
            break;

        case "vcf_consensus":
            check_args(consensus,
                       consensus_args, 
                       vcf_consensus_accept_args)
            break;

        default: 
           check_args(consensus, 
                      consensus_args,
                      bcftools_accept_args)
    } 
}
//----------------------------------------

//----------------------------------------
// HELPER FUNCTIONS: CHECK FILTER ARGS 
//----------------------------------------
def check_filter(filter) {
    println("Checking filter tool arguments")
    filter_args = [:]
    arg_str_to_dict(params.filter_args, filter_args)
    
    bcftools_accept_args = ['key', 'bcftools']
    
     switch(filter) {
        case "bcftools": 
           check_args(filter,
                      filter_args,
                      bcftools_accept_args)
            break;

        default: 
            check_args(filter,
                       filter_args,
                       bcftools_accept_args)
    }
}
//----------------------------------------

//=============================================================================
// INPUT VALIDATION
//=============================================================================

input_validation('align', params.align, 'align_args')
input_validation('variant', params.variant, 'variant_args')
input_validation('filter', params.filter, 'filter_args')
input_validation('consensus', params.consensus, 'consensus_args')

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
summary['Filter']               = params.filter
summary['Filter Args']          = params.filter_args

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

//----------------------------------------
// PROCESSES: RESULTDIR
//----------------------------------------
process RESULTDIR {
    publishDir "${params.outdir}/align/alignment", 
                pattern: "*.txt", mode: 'copy'
     
    output: 
        file "readme.txt"
    
    script:
        alignment = "readme.txt"
    
    """
    echo "readme" > $alignment
    """
}
//----------------------------------------

//----------------------------------------
// PROCESSES: REMOVE_PHIX
// Remove any Coliphage phi-X714 reads using bbduk
// bbduk and options explained above script call
//----------------------------------------
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
//----------------------------------------

//----------------------------------------
// PROCESSES: FASTP
// Adapter trimming and quality filering using fastp
// By default, only up to 40% of bases can be below the min Phred score base
// quality threshold of Q15
// fastp and options explained above script call
//----------------------------------------
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
//----------------------------------------

//----------------------------------------
// PROCESSES: FASTQC
// Quality control using fastqc
// fastqc and options explained above script call
//----------------------------------------
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
//----------------------------------------

//----------------------------------------
// PROCESSES: UNZIP 
// Unzip file pair 
//----------------------------------------

process UNZIP {
    tag "$sample_id"
    publishDir "${params.outdir}/align/bwa",
                pattern: "*.fastp.fastq", mode: 'copy'
    
    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path(r1_out), path(r2_out),
                emit: 'reads'

    script:
        r1_out = "${sample_id}_1.fastp.fastq"
        r2_out = "${sample_id}_2.fastp.fastq"

    """
    gzip -d --force $r1 > $r1_out;
    gzip -d --force $r2 > $r2_out;
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: BWAINDEX 
// Build the reference genome index
//----------------------------------------

process BWAINDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/align/bwa", 
                pattern: "*.sai", mode: 'copy'
    publishDir "${params.outdir}/align/bwa/index",
                pattern: "ref.fa.*", mode: 'copy'
    input:
        file ref            
        tuple val(sample_id), path(r1), path(r2)
    
    output:
        file "${sample_id}_1.sai"
        file "${sample_id}_2.sai"
        tuple val(sample_id), path(r1), path(r2), path(r1_sai), path(r2_sai),
                emit: 'reads'
    
    script:
        r1_sai = "${sample_id}_1.sai"
        r2_sai = "${sample_id}_2.sai"

    """
    bwa index $ref;
    bwa aln $ref $r1 > $r1_sai;
    bwa aln $ref $r2 > $r2_sai;
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: BWA 
// Paired End read alignment 
// Utilizes reference genome
//----------------------------------------

process BWA {
    tag "$sample_id"
    publishDir "${params.outdir}/align/alignment",
                pattern: "*.sam", mode: 'copy'

    input:
        file ref            
        tuple val(sample_id), path(r1), path(r2), path(r1_sai), path(r2_sai)
    
    output:
        file "${sample_id}_bwa_align_pe.sam"
        val sample_id, emit: 'sample_id'
    
    script:
        align = "${sample_id}_bwa_align_pe.sam"
    
    """
    bwa index $ref;
    bwa sampe $ref $r2_sai $r2_sai $r1 $r2 > $align;
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: BOWTIE2 
// Build the reference genome index
//----------------------------------------

process BOWTIE2 {
    tag "$sample_id"
    publishDir "./index", pattern: "index*bt2*", mode: "copy"
    publishDir "${params.outdir}/align/alignment",
                pattern: "*.sam", mode: 'copy'
    input:
        file ref
        tuple val(sample_id), path(r1), path(r2)

    output:
        file "${sample_id}_bowtie2_align_pe.sam"
        val sample_id, emit: 'sample_id'

    script:
        index_dir = "./index"
        indexes = "./index/index"
        align = "${sample_id}_bowtie2_align_pe.sam"
    
    """
    sudo mkdir $index_dir;
    sudo chmod 777 $index_dir;
    bowtie2-build $ref $indexes;
    bowtie2 -x $indexes -1 $r1 -2 $r2 > $align;
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: FASTQTOFASTA
// Convert fastq files to fasta files 
//----------------------------------------

process FASTQTOFASTA {
    tag "$sample_id"
    publishDir "${params.outdir}/align/minimap2",
                pattern: "*.fasta", mode: 'copy'

    input:
        tuple val(sample_id), path(r1), path(r2)
    
    output:
        tuple val(sample_id), path(r1_fasta), path(r2_fasta),
                emit: 'reads'
    
    script:
        r1_fasta = "${sample_id}_1.fastp.fasta"
        r2_fasta = "${sample_id}_2.fastp.fasta"
    
    """
    sed -n '1~4s/^@/>/p;2~4p' $r1 > $r1_fasta;
    sed -n '1~4s/^@/>/p;2~4p' $r2 > $r2_fasta;

    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: MINIMAP2 
// Paired End read alignment 
// Utilizes reference genome
//----------------------------------------

process MINIMAP2 {
    tag "$sample_id"
    publishDir "${params.outdir}/align/alignment",
                pattern: "*.sam", mode: 'copy'

    input:
        file ref            
        tuple val(sample_id), path(r1), path(r2)
    
    output:
        file "${sample_id}_minimap2_align_pe.sam"
        val sample_id, emit: 'sample_id'
    
    script:
        align = "${sample_id}_minimap2_align_pe.sam"

    """
    minimap2 -ax sr $ref $r1 $r2 > $align;
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: SAMTOBAM
// Convert files from SAM to BAM format  
// Prepares output from BWA to become input for Freebayes
//----------------------------------------

process SAMTOBAM {
    tag "$sample_id"
    publishDir "${params.outdir}/samconvert/${method}",
                pattern: "*.bam", mode: "copy"

    input: 
        file align
        val sample_id

    output:
        tuple val(sample_id), path(align_bam), val(method), 
                emit: 'align'
        
    script: 
        method = "${align}".split('_')[1]
        align_bam = "${sample_id}_${method}_align_pe.bam"

    """
    samtools view -S -b $align > $align_bam;
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: BAMSORT
// Sort Bam file  
//----------------------------------------

process BAMSORT {
    tag "$sample_id"
    publishDir "${params.outdir}/samconvert/${method}",
                pattern: "*.sorted.bam", mode: "copy"

    input: 
        tuple val(sample_id), path(align), val(method)

    output:
        tuple val(sample_id), path(sorted_bam), 
                emit: 'align'
    
    script: 
        sorted_bam = "${sample_id}_${method}_align_pe.sorted.bam"

    """
    samtools sort $align -o $sorted_bam -O bam
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: BAMINDEX
// Index sorted Bam file  
//----------------------------------------

process BAMINDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/samconvert",
                pattern: "*.bai", mode: "copy"

    input: 
        tuple val(sample_id), path(align)

    output:
        tuple val(sample_id),  path(align),
                emit: 'align'
    
    script: 
        indexed_bam = "${sample_id}_align_pe.bai"

    """
    samtools index $align $indexed_bam
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: FREEBAYES
// Variant Calling 
// Freebayes is the default variant caller  
//----------------------------------------

process FREEBAYES {
    tag "$sample_id"
    publishDir "${params.outdir}/variant/freebayes",
                pattern: "*.vcf", mode: 'copy'

    input: 
        file ref
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path(variant), 
                emit: 'variant'
    
    script: 
        variant = "${sample_id}_var.vcf"

    """
    freebayes -f $ref $bam > $variant
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: BCFTOOLS_STATS
// Collect stats from variant calls 
//----------------------------------------

process BCFTOOLS_STATS {
    tag "$sample_id"
    publishDir "${params.outdir}/variant",
                pattern: "*_stats.txt", mode: 'copy'

    input: 
        tuple val(sample_id), path(variant)

    output:
        file "${sample_id}_stats.txt"
            
    script: 
        stats = "${sample_id}_stats.txt"

    """
    bcftools stats $variant > $stats; 
    
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: BCFTOOLS_FILTER
// Filtering variant calls 
//----------------------------------------

process BCFTOOLS_FILTER {
    tag "$sample_id"
    publishDir "${params.outdir}/variant",
                pattern: "*_filtered.vcf", mode: 'copy'

    input: 
        tuple val(sample_id), path(variant)

    output:
        tuple val(sample_id), path(variant),
                emit: 'variant'
        file "${sample_id}_filtered.vcf"
            
    script: 
        filtered = "${sample_id}_filtered.vcf"

    """
    bcftools stats $variant > $filtered; 
    
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: SNPEFF
// Annotates variants 
//----------------------------------------

process SNPEFF {
    tag "$sample_id"
    publishDir "${params.outdir}/variant/prediction",
                pattern: "*_annotation.vcf", mode: 'copy'

    input: 
        tuple val(sample_id), path(variant)

    output:
        file "${sample_id}_annotation.vcf"
         tuple val(sample_id), path(variant), 
                emit: 'annotation' 
    
    script:
        filtered = "${sample_id}_annotation.vcf"

    """
    java -jar /usr/local/src/snpEff/snpEff.jar -v 
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: BCFTOOLS_CONSENSUS
// Build consensus sequence from variants
//----------------------------------------

process BCFTOOLS_CONSENSUS {
    tag "$sample_id"
    publishDir "${params.outdir}/consensus/bcftools",
                pattern: "*.txt", mode: 'copy'

    input: 
        tuple val(sample_id), path(variant)

    output:
       
    script:

    """
    bcftools consensus
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: VCF_CONSENSUS
// Build consensus sequence from variants
//----------------------------------------

process VCF_CONSENSUS {
    tag "$sample_id"
    publishDir "${params.outdir}/consensus/vcf",
                pattern: "*.txt", mode: 'copy'

    input:
        tuple val(sample_id), path(variant)

    output:

    script:

    """
    bcftools consensus
    """
}

//----------------------------------------

//=============================================================================
// WORKFLOW DEFINITION 
//=============================================================================

workflow {
    //----------------------------------------
    // Read input, trimming and quality control
    //----------------------------------------
    // Create channel for paired end reads
    println("workflow") 
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
    
    
    //----------------------------------------
    // Aligners Put each aligner in a different workflow
    //----------------------------------------
    ch_ref = Channel.value(file("${baseDir}/data/ref.fa"))
    ch_align = Channel.watchPath("${params.outdir}/align/alignment/*.sam")
                        .view{ "$it" }

    UNZIP(FASTP.out.reads)
    
    if (params.align == 'bwa') {     
        
        BWAINDEX(ch_ref, UNZIP.out.reads)
        BWA(ch_ref, BWAINDEX.out.reads)
        SAMTOBAM(ch_align, BWA.out.sample_id)
    
    } else if (params.align == 'bowtie2') {
        
        BOWTIE2(ch_ref, UNZIP.out.reads)
        SAMTOBAM(ch_align, BOWTIE2.out.sample_id)

    } else if (params.align == 'minimap2') {
        
        FASTQTOFASTA(UNZIP.out.reads)
        MINIMAP2(ch_ref, FASTQTOFASTA.out.reads)
        SAMTOBAM(ch_align, MINIMAP2.out.sample_id)
    }
    
    //----------------------------------------
    // From samtools to SnpEff
    //----------------------------------------
    BAMSORT(SAMTOBAM.out.align)
    BAMINDEX(BAMSORT.out.align)

    //freebayes process
    FREEBAYES(ch_ref, BAMINDEX.out.align)
     
    BCFTOOLS_STATS(FREEBAYES.out.variant)
    BCFTOOLS_FILTER(FREEBAYES.out.variant)
    SNPEFF(BCFTOOLS_FILTER.out.variant)
    
    //----------------------------------------
    // Consensus building put in own workflows like aligners
    //----------------------------------------



}

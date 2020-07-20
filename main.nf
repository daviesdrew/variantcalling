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
summary['Reference']            = params.ref

summary['Aligner']              = params.align
summary['Aligner Args']         = params.align_args
summary['Variant Caller']       = params.variant
summary['Variant Caller Args']  = params.variant_args
summary['Filter']               = params.filter
summary['Filter Args']          = params.filter_args
summary['Consensus']            = params.consensus
summary['Consensus Args']       = params.consensus_args

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

log.info summary.collect { k, v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


//=============================================================================
// WORKFLOW DEFINITION 
//=============================================================================
include bwa from "./modules/pipes.nf"
include bowtie2 from "./modules/pipes.nf"
include minimap2 from "./modules/pipes.nf"

include quality_check from "./modules/quality_check.nf"

workflow {
    main:    
        ref = Channel.fromPath(params.ref) 
        phix = Channel.value(file("${baseDir}/${params.phix}"))
        reads = Channel.fromFilePairs(params.reads,
                                      flat: true,
                                      checkIfExists: true)
                       //check if Channel is empty
                       .ifEmpty{ exit 1, 
                                 """No reads specified! 
                                    Please specify where your reads are, 
                                    
                                    e.g. '--reads \"/path/to/reads/*R{1,2}.fastq.qz\"' 
                            
                                    (quotes around reads ath required if using `*` and 
                                     other characters expanded by the shell!)"""}
                       //Convert specified files into tuples
                       .map{ [ it[0].replaceAll(/_S\d{1,2}_L001/,""), it[1], it[2] ] }
        
        quality_check(reads, phix)
        if(params.align == 'bowtie2' || params.align == 'all')
            bowtie2(ref, quality_check.out.reads)

        if(params.align == 'minimap2' || params.align == 'all')
            minimap2(ref, quality_check.out.reads) 
        
        if(params.align == 'bwa' || params.align == 'all')
            bwa(ref, quality_check.out.reads)

}

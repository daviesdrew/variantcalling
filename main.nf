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

//check for installation of bwa aligner
//ensure arguments are valid for bwa
def check_bwa(bwa, args) {

}

//check for installation of minimap2 aligner
//ensure arguments are valid for minimap2
def check_minimap2(minimap2, args) {

}

//check for installation of bowtie2 aligner
//ensure arguments are valid for bowtie2
def check_bowtie2(bowtie2, args) {

}

//check for installation of aligner
//check for argument validitiy
def check_aligner(aligner, align_args) {
    switch(aligner) {
        case "bwa": 
            println("BWA aligner selected");
            check_bwa(aligner, align_args);
            break;

        case "minimap2":
            println("Minimap2 aligner selected");
            check_minimap2(aligner, align_args);
            break;

        case "bowtie2": 
            println("Bowtie2 aligner selected");
            check_bowtie2(aligner, align_args);
            break;

        default: 
            println("Default aligner selected");
            println("Default aligner is set to Bowtie2");
            check_bowtie2(aligner, align_args);
    }
}

//=============================================================================
// INPUT VALIDATION
//=============================================================================

align_args = params.align_args.tokenize(',')

align_args.args.each{ k -> println "${k}" }

println "${align_args}"

check_aligner(params.aligner, params.align_args)


//=============================================================================
// WORKFLOW RUN PARAMETERS LOGGING
//=============================================================================

//=============================================================================
// PROCESSES
//=============================================================================

//=============================================================================
// WORKFLOW DEFINITION 
//=============================================================================



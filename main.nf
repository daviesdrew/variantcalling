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
    
    
    
    
    ${c_bul}BWA Aligner Options:${c_reset}
      --max_edit_dist []             Desc ': 'n', //bwa aln -n
                    
      --max_gap_opens []             Desc ': 'o', //bwa aln -o
                    
      --max_gap_ext []               Desc ': 'e', //bwa aln -e
      
      --no_long_del []               Desc ': 'd', //bwa aln -d 
      
      --limit_indel []               Desc ': 'i', //bwa aln -i
      
      --subseq_seed []               Desc ': 'l', //bwa aln -l
      
      --mismatch_pen []              Desc ': 'M', //bwa aln -M
      
      --gap_open_pen []              Desc ': 'O', //bwa aln -O
      
      --max_seed_edit_dist []        Desc ': 'k', //bwa aln -k
      
      --gap_ext_pen []               Desc ': 'E', //bwa aln -E
      
      --max_insert []                Desc ': 'a', //bwa sampe -a
      
      --max_occur []                 Desc ': 'o', //bwa sampe -o
      
      --max_align []                 Desc ': 'n', //bwa sampe -n
      
      --max_discord []               Desc ': 'N', //bwa sampe -N
      
      --index_algo []                Desc ': 'a' //bwa index -a

   ${c_bul}Bowtie2 Aligner Options:${c_reset}
     --end_to_end []                 Desc 
     
     --very_fast_local []            Desc
     
     --fast_local []                 Desc
     
     --sensitive_local []            Desc 
    
     --very_sensitive_local []       Desc

     --dpad []                       Desc

     --gbar []                       Desc

     --ext_attempt_count []          Desc

     --max_reseed []                 Desc

     --minins []                     Desc

     --maxins []                     Desc

     --preseq_assay []               Desc

     --no_mixed []                   Desc

     --no_contain []                 Desc

     --no_lap []                     Desc

     --qc_filter []                  Desc

   ${c_bul}Minimap2 Aligner Options:${c_reset}
     --minr_kmer_len []              Desc 
     
     --minr_win_size []              Desc
     
     --homopoly_min []               Desc
     
     --save_index []                 Desc 
    
     --ignore_minr []                Desc

     --min_occ_floor []              Desc

     --stop_chain_elong []           Desc

     --approx_gap_size []            Desc

     --discard_chain_minr []         Desc

     --discard_chain_score []        Desc

     --max_ref_gap []                Desc

     --max_frag_len []               Desc

     --splice []                     Desc

     --sr []                         Desc

     --for_only []                   Desc

     --rev_only []                   Desc

     --no_mixed []                   Desc

     --no_contain []                 Desc

     --no_lap []                     Desc

     --qc_filter []                  Desc

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

//=============================================================================
// CONSTANT VALUES
//=============================================================================

//----------------------------------------
// CONSTANT VALUES: ALIGNERS
//----------------------------------------

params.bwa_args = [ 
    'max_edit_dist': 'n', //bwa aln -n
    'max_gap_opens': 'o', //bwa aln -o
    'max_gap_ext': 'e', //bwa aln -e
    'no_long_del': 'd', //bwa aln -d 
    'limit_indel': 'i', //bwa aln -i
    'subseq_seed': 'l', //bwa aln -l
    'mismatch_pen': 'M', //bwa aln -M
    'gap_open_pen': 'O', //bwa aln -O
    'max_seed_edit_dist': 'k', //bwa aln -k
    'gap_ext_pen': 'E', //bwa aln -E
    'max_insert': 'a', //bwa sampe -a
    'max_occur': 'o', //bwa sampe -o
    'max_align': 'n', //bwa sampe -n
    'max_discord': 'N', //bwa sampe -N
    'index_algo': 'a' //bwa index -a
]

params.bowtie2_args = [ 
    'end_to_end': '-end-to-end', //bowtie2 --end-to-end
    'very_fast_local': '-very-fast-local', //bowtie2 --very-fast-local
    'fast_local': '-fast-local', //bowtie2 --fast-local
    'sensitive_local': '-sensitive-local', //bowtie2 --sensitive-local
    'very_sensitive_local': '-very-sensitive-local', //bowtie2 --very-sensitive-local
    'dpad': '-dpad', //bowtie2 --dpad
    'gbar': '-bgar', //bowtie2 --gbar
    'ext_attempt_count': 'D', //bowtie2 -D [int]
    'max_reseed': 'R', //bowtie2 -R [int]
    'minins': '-minins', //bowtie2 --minins [int]
    'maxins': '-maxins', //bowtie2 --maxins [int]
    'preseq_assay': '-fr', //bowtie2 --fr
    'no_mixed': '-no-mixed', //bowtie2 --no-mixed
    'no_contain': '-no-contain', //bowtie2 --no-contain
    'no_lap': '-no-overlap', //bowtie2 --no-overlap
    'qc_filter': '-qc_filter', //bowtie2 --qc-filter
]

params.minimap2_args = [ 
    'minr_kmer_len': 'k', //minimap2 -k [int]
    'minr_win_size': 'w', //minimap2 -w [int]
    'homopoly_minr': 'H', //minimap2 -H
    'save_index': 'd', //minimap2 -d [file]
    'ignore_minr': 'f', //minimap2 -f [float|int]
    'min_occ_floor': '-min-occ-floor', //minimap2 --min-occ-floor [int]
    'stop_chain_elong': 'g', //minimap2 -g [int]
    'approx_gap_size': 'r', //minimap2 -r [int]
    'discard_chain_minr': 'n', //minimap2 -n [int]
    'discard_chain_score': 'm', //minimap2 -m [int]
    'num_sec_align': 'N', //minimap2 -N [int]
    'max_ref_gap': 'G', //minimap2 -G [num]
    'max_frag_len': 'F', //minimap2 -F [num]
    'splice': '-splice', //minimap2 --splice
    'sr': '-sr', //minimap2 --sr
    'for_only': '-for-only', //minimap2 --for-only
    'rev_only': '-rev-only', //minimap2 --rev-only
    'match_score': 'A', //minimap2 -A [int]
    'mismatch_pen': 'B', //minimap2 -B [int]
    'gap_open_pen': 'O', //minimap2 -O [int]
    'gap_ext_pen': 'E', //minimap2 -E [int]
    'non_canonical': 'C', //minimap2 -C [int]
    'end_bonus': '-end-bonus', //minimap2 --end-bonus [int]
    'ambig_mismatch': '-score-N' //minimap2 --score-N [int]
]

//----------------------------------------

//----------------------------------------
// CONSTANT VALUES: VARIANT CALLERS
//----------------------------------------

params.freebayes_args = [ 
    'key': 'one' 
]

//----------------------------------------

//----------------------------------------
// CONSTANT VALUES: FILTERS
//----------------------------------------

params.bcftools_filter_args = [
    'key': 'one'
]
//----------------------------------------

//----------------------------------------
// CONSTANT VALUES: CONSENSUSES 
//----------------------------------------

params.bcftools_consensus_args = [
    'key': 'one'
]

params.vcf_consensus_args = [
    'key': 'one'
]

//----------------------------------------

//=============================================================================
// HELPER FUNCTIONS
//=============================================================================

//----------------------------------------
// HELPER FUNCTIONS: GENERAL TOOLS
//----------------------------------------
def check_arg_existence(tool, tool_args) {
    params[tool] = (params.containsKey(tool)) ?
                        params[tool] :
                        null
    params[tool_args] = (params.containsKey(tool_args)) ?
                        params[tool_args] :
                        ""
}

def check_dict_args(args, acceptable_args) {
    args.keySet()
        .each{ arg -> assert acceptable_args.contains(arg) }
    
    return args
}

def swap_arg_keys(pipeline_args, cmd_args, new_dict) {
    pipeline_args.keySet()
                 .each{ key -> new_dict[cmd_args[key]] = pipeline_args[key] }
    return new_dict 
}

def arg_str_to_dict(arg_str, dict) {
    arg_str.tokenize(',')
        .each{ arg -> 
            dict[arg.tokenize()[0].replace('-','')] = arg.tokenize()[1] }
}

def arg_dict_to_str(dict) {
    str = ""
    dict.each{ k, v -> str = str.concat('-'+k+' ').concat(v+' ') }
    return str
}

def print_arg_comparison(tool, args, acceptable_args) {
    println("${tool} selected\n")
    println("${tool} args")
    println("${args}\n")
    println("${tool} accepted args")
    println("${acceptable_args}\n")
}

def check_args(tool, tool_args, tool_accept_args) {
    tool_arg_keys = tool_accept_args.keySet()
    print_arg_comparison(tool, tool_args, tool_arg_keys)
    args = check_dict_args(tool_args, tool_arg_keys)
    new_dict = [:]
    new_dict = swap_arg_keys(tool_args, tool_accept_args, new_dict)
    //println("new_dict")
    //println("${new_dict}")
    //println("new string from dict")
    str = arg_dict_to_str(new_dict)
    //println("${str}")   
    return str
}


//----------------------------------------
// HELPER FUNCTIONS: INPUT VALIDATION
//----------------------------------------

def input_validation(tool_type, tool_args) {
    tool = params[tool_type]
    args = null
    switch(tool_type) {
            
            case 'align':
                args = check_aligner(tool)
                //check_sys_for_aligner(tool)
                break;

            case 'variant': 
                args = check_variant_caller(tool)
                //check_sys_for_variant_caller(tool)
                break;

            case 'filter':
                args = check_filter(tool)
                //check_sys_for_filter(tool)
                break;
        
            case 'consensus':
                args = check_consensus(tool)
                //check_sys_for_consensus(tool)
                break;
    }
    println("Args with cmd specific options")
    println("${args}")

}
//----------------------------------------

//----------------------------------------
// HELPER FUNCTIONS: CHECK ALIGNER ARGS
//----------------------------------------
def check_aligner(aligner) {
    println("Checking aligner arguments")

    align_args = [:]
    arg_str_to_dict(params.align_args, align_args)

    
    acceptable_args = [ 
        'bwa': params.bwa_args,
        'minimap2': params.minimap2_args,
        'bowtie2': params.bowtie2_args     
    ]

    if (aligner == null) 
        aligner = 'bwa'

    check_args(aligner, 
               align_args,
               acceptable_args[aligner])
    return str
}
//----------------------------------------

//----------------------------------------
// HELPER FUNCTIONS: CHECK VARIANT CALLING ARGS
//----------------------------------------
def check_variant_caller(variant_caller) {
    println("Checking variant caller arguments")
    variant_caller_args = [:]
    arg_str_to_dict(params.variant_args, variant_caller_args)
    
    acceptable_args = [
        'freebayes': params.freebayes_args
    ]
    
    if(variant_caller == null) 
        variant_valler = 'freebayes'

    return check_args(variant_caller,
                      variant_caller_args,
                      acceptable_args[variant_caller])
}
//----------------------------------------

//----------------------------------------
// HELPER FUNCTIONS: CHECK FILTER ARGS 
//----------------------------------------
def check_filter(filter) {
    println("Checking filter tool arguments")
    filter_args = [:]
    arg_str_to_dict(params.filter_args, filter_args)
    
    acceptable_args = [
        'bcftools': params.bcftools_filter_args
    ]
    
    return check_args(filter,
                      filter_args,
                      acceptable_args[filter])
}

//----------------------------------------

//----------------------------------------
// HELPER FUNCTIONS: CHECK CONSENSUS ARGS
//----------------------------------------
def check_consensus(consensus) {
    println("Checking consensus arguments")
    consensus_args = [:]
    
    if (consensus != null) {
        arg_str_to_dict(params.consensus_args, consensus_args)
    }

    acceptable_args = [
        'bcftools': params.bcftools_consensus_args,
        'vcf_consensus': params.vcf_consensus_args
    ]

    if (consensus == null)
        consensus = 'bcftools'
    
    return check_args(consensus,
                      consensus_args,
                      acceptable_args[consensus])
}
//----------------------------------------

//=============================================================================
// INPUT VALIDATION
//=============================================================================

arguments = [ 'align': 'align_args', 
              'variant': 'variant_args',
              'filter': 'filter_args',
              'consensus': 'consensus_args' ] 

arguments.each{ k, v -> check_arg_existence(k, v) }
arguments.each{ k, v -> input_validation(k, v) }

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
summary['Consensus']            = (params.containsKey('consensus'))
                                    ? params.consensus
                                    : null
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

//=============================================================================
// PROCESSES
//=============================================================================

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
    bwa index -a bwtsw $ref;
    bwa aln -t ${task.cpus} -I $ref $r1 > $r1_sai;
    bwa aln -t ${task.cpus} -I $ref $r2 > $r2_sai;
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
        tuple val(sample_id), path(align), emit: 'align'
    
    echo true
    
    script:
        align = "${sample_id}_bwa_align_pe.sam"
        args = "${params.align_args}"
      
    """
    echo $args
    bwa index -a bwtsw $ref;
    bwa sampe -P $ref $r2_sai $r2_sai $r1 $r2 > $align;
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
        tuple val(sample_id), path(align), emit: 'align'

    script:
        index_dir = "./index"
        indexes = "./index/index"
        align = "${sample_id}_bowtie2_align_pe.sam"
    
    """
    sudo mkdir $index_dir;
    sudo chmod 777 $index_dir;
    bowtie2-build $ref $indexes;
    bowtie2 --threads ${task.cpus} -x $indexes -1 $r1 -2 $r2 -S $align;
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
        tuple val(sample_id), path(align), emit: 'align'
    
    script:
        align = "${sample_id}_minimap2_align_pe.sam"

    """
    minimap2 -a -t ${task.cpus} -ax sr $ref $r1 $r2 -o $align;
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
        tuple val(sample_id), path(align) 

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
        tuple val(sample_id), path(sorted_bam), val(method),
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
        tuple val(sample_id), path(align), val(method)

    output:
        tuple val(sample_id),  path(align), val(method),
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
        tuple val(sample_id), path(bam), val(method)

    output:
        tuple val(sample_id), path(variant), path(ref), 
                emit: 'variant'
    
    script: 
        variant = "${sample_id}.vcf"

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
        tuple val(sample_id), path(variant), path(ref)

    output:
        file "${sample_id}_stats.txt"
            
    script: 
        stats = "${sample_id}_stats.txt"

    """
    bcftools stats -F $ref $variant > $stats; 
    
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
        tuple val(sample_id), path(variant), path(ref)

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
        annotated = "${sample_id}_annotation.vcf"

    """
    java -jar /usr/local/src/snpEff/snpEff.jar asfv $variant -v \\
    > $annotated
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
                pattern: "*.vcf", mode: 'copy'

    input: 
        file ref 
        tuple val(sample_id), path(variant)

    output:
        file "${sample_id}_consensus.vcf"
       
    script:
        zipped = "${variant}.gz"
        consensus = "${sample_id}_consensus.vcf"

    """
    bgzip $variant
    tabix -p vcf $zipped
    bcftools consensus -f $ref $zipped > $consensus
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

    UNZIP(FASTP.out.reads)
    
    if (params.align == 'bwa') {     
        
        BWAINDEX(ch_ref, UNZIP.out.reads)
        BWA(ch_ref, BWAINDEX.out.reads)
        SAMTOBAM(BWA.out.align) //ch_align, BWA.out.sample_id)
    
    } else if (params.align == 'bowtie2') {
        
        BOWTIE2(ch_ref, UNZIP.out.reads)
        SAMTOBAM(BOWTIE2.out.align)

    } else if (params.align == 'minimap2') {
        
        FASTQTOFASTA(UNZIP.out.reads)
        MINIMAP2(ch_ref, FASTQTOFASTA.out.reads)
        SAMTOBAM(MINIMAP2.out.align)
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
    
    if (params.consensus == 'bcftools') {

        BCFTOOLS_CONSENSUS(ch_ref, SNPEFF.out.annotation)

    } else if (params.consensus == 'vcf_consensus') {

        VCF_CONSENSUS(ch_ref, SNPEFF.out.annotation)
    
    } else {

        BCFTOOLS_CONSENSUS(ch_ref, SNPEFF.out.annotation)
    
    }
}

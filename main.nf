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
      --max_edit_dist []             Desc bwa aln -n
                    
      --max_gap_opens []             Desc bwa aln -o
                    
      --max_gap_ext []               Desc bwa aln -e
      
      --no_long_del []               Desc bwa aln -d 
      
      --limit_indel []               Desc bwa aln -i
      
      --subseq_seed []               Desc bwa aln -l
      
      --mismatch_pen []              Desc bwa aln -M
      
      --gap_open_pen []              Desc bwa aln -O
      
      --max_seed_edit_dist []        Desc bwa aln -k
      
      --gap_ext_pen []               Desc bwa aln -E
     
      --max_insert []                Desc bwa sampe -a
      
      --max_occur []                 Desc bwa sampe -o
      
      --max_align []                 Desc bwa sampe -n
      
      --max_discord []               Desc bwa sampe -N
      
      --index_algo []                Desc bwa index -a

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

     --match_score []                Desc

     --mismatch_pen []               Desc

     --gap_open_pen []               Desc

     --gap_ext_pen []                Desc

     --non_canonical []              Desc

     --end_bonus []                  Desc

     --ambig_mismatch []             Desc


   ${c_bul}Freebayes Variant Caller Options:${c_reset}
     --use_ref_allele []             Desc freebayes --use-reference-allele
     
     --theta []                      Desc freebayes --theta [num]
    
     --ploidy []                     Desc freebayes --ploidy [int]
    
     --nbest_allele []               Desc freebayes -use-best-n-alleles [int]
     
     --max_cmplx_gap []              Desc freebayes --max-complex-gap [int]
    
     --haplo_len []                  Desc freebayes --haplotype-length [int]

     --min_rep_size []               Desc freebayes --min-repeat-size [int]
    
     --min_rep_entropy []            Desc freebayes --min-repeat-entropy [int]
    
     --no_part_obs []                Desc freebayes --no-partial-observations 
    
     --use_dup_reads []              Desc freebayes --use-duplicate-reads
    
     --base_qual_cap []              Desc freebayes --base-quality-cap [Q]
    
     --min_map_qual []               Desc freebayes --min-mapping-quality [Q]
    
     --min_base_qual []              Desc freebayes --min-base-quality [Q]
    
     --mismatch_qual_thresh []       Desc freebayes --mismatch-base-quality-threshold [Q]
    
     --mismatch_limit []             Desc freebayes --read-mismatch-limit [int]
    
     --read_snp_limit []             Desc freebayes --read-snp-limit [int]
    
     --read_indel_limit []           Desc freebayes --read-indel-limit [int]

   ${c_bul}Freebayes Variant Caller Options:${c_reset}
     --soft_filter                   Desc bcftools filter --soft-filter [string]
    
    --set_GTs                        Desc bcftools filter --set-GTs [.|0]
    
    --snp_gap                        Desc bcftools filter --SnpGap [int]

    --indel_gap                      Desc bcftools filter --IndelGap [int]

    --mode                           Desc bcftools filter --mode [+x]

   ${c_bul}SnpEff Gene Prediction Options:${c_reset}
    
    
    
    
    
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
    'use_ref_allele': '-use-reference-allele', //freebayes --use-reference-allele
    'theta': '-theta', //freebayes --theta [num]
    'ploidy': '-ploidy', //freebayes --ploidy [int]
    'nbest_allele': '-use-best-n-alleles', //freebayes -use-best-n-alleles [int]
    'max_cmplx_gap': '-mask-complex-gap', //freebayes --max-complex-gap [int]
    'haplo_len': '-haplotype-length', //freebayes --haplotype-length [int]
    'min_rep_size': '-min-repeat-size', //freebayes --min-repeat-size [int]
    'min_rep_entropy': '-min-repeat-entropy', //freebayes --min-repeat-entropy [int]
    'no_part_obs': '-no-partial-observations', //freebayes --no-partial-observations 
    'use_dup_reads': '-use-duplicate-reads', //freebayes --use-duplicate-reads
    'base_qual_cap': '-base-quality-cap', //freebayes --base-quality-cap [Q]
    'min_map_qual': '-min-mapping-quality', //freebayes --min-mapping-quality [Q]
    'min_base_qual': '-min-base-quality', //freebayes --min-base-quality [Q]
    'mismatch_qual_thresh': '-mismatch-base-quality-threshold', //freebayes --mismatch-base-quality-threshold [Q]
    'mismatch_limit': '-read-mismatch-limit', //freebayes --read-mismatch-limit [int]
    'read_snp_limit': '-read-snp-limit', //freebayes --read-snp-limit [int]
    'read_indel_limit': '-read-indel-limit' //freebayes --read-indel-limit [int]
]

//----------------------------------------

//----------------------------------------
// CONSTANT VALUES: FILTERS
//----------------------------------------

params.bcftools_filter_args = [
    'include': '-include', //bcftools filter --include [expr]

    'exclude': '-exclude', //bcftools filter --exclude [expr]

    'soft_filter': '-soft-filter', //bcftools filter --soft-filter [string]
    
    'set_GTs': '-set-GTs', //bcftools filter --set-GTs [.|0]
    
    'snp_gap': '-SnpGap', //bcftools filter --SnpGap [int]

    'indel_gap': '-IndelGap', //bcftools filter --IndelGap [int]

    'mode': '-mode' //bcftools filter --mode [+x]
]
//----------------------------------------

//----------------------------------------
// CONSTANT VALUES: PREDICTION
//----------------------------------------

params.snpeff_args = [
    'soft_filter': '-soft-filter', // bcftools filter --soft-filter [string]
    
    'set_GTs': '-set-GTs', //bcftools filter --set-GTs [.|0]
    
    'snp_gap': '-SnpGap', //bcftools filter --SnpGap [int]

    'indel_gap': '-IndelGap', //bcftools filter --IndelGap [int]

    'mode': '-mode' //bcftools filter --mode [+x]
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
    return arg_dict_to_str(new_dict)
}


//----------------------------------------
// HELPER FUNCTIONS: INPUT VALIDATION
//----------------------------------------

def input_validation(tool_type, tool_args) {
    tool = params[tool_type]
    args = null
    switch(tool_type) {
            case 'align':
                check_aligner(tool)
                println("Args for align process")
                println("${params.align_args_str}\n")
                //check_sys_for_aligner(tool)
                break;

            case 'variant': 
                check_variant_caller(tool)
                println("Args for variant calling process")
                println("${params.variant_args_str}\n")
                //check_sys_for_variant_caller(tool)
                break;

            case 'filter':
                check_filter(tool)
                println("Args for filter process")
                println("${params.filter_args_str}\n")
                //check_sys_for_filter(tool)
                break;
        
            case 'prediction':
                check_prediction(tool)
                println("Args for prediction process")
                println("${params.prediction_args_str}\n")
                //check_sys_for_prediction(tool)
                break;

            case 'consensus':
                check_consensus(tool)
                println("Args for consensus process")
                println("${params.consensus_args_str}\n")
                //check_sys_for_consensus(tool)
                break;
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

    
    acceptable_args = [ 
        'bwa': params.bwa_args,
        'minimap2': params.minimap2_args,
        'bowtie2': params.bowtie2_args     
    ]

    if (aligner == null) 
        aligner = 'bwa'

    params.align_args_str = check_args(aligner, 
                                       align_args,
                                       acceptable_args[aligner])
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

    params.variant_args_str = check_args(variant_caller,
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
    
    params.filter_args_str = check_args(filter,
                                        filter_args,
                                        acceptable_args[filter])
}

//----------------------------------------

//----------------------------------------
// HELPER FUNCTIONS: CHECK PREDICTION ARGS 
//----------------------------------------
def check_prediction(prediction) {
    println("Checking prediction tool arguments")
    prediction_args = [:]
    arg_str_to_dict(params.prediction_args, prediction_args)
    
    acceptable_args = [
        'snpeff': params.snpeff_args
    ]
    println("before check_args") 
    print_arg_comparison(prediction, 
                          prediction_args, 
                          acceptable_args[prediction])
     
    params.prediction_args_str = check_args(prediction,
                                            prediction_args,
                                            acceptable_args[prediction])
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
    
    params.consensus_args_str = check_args(consensus,
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
              'prediction': 'prediction_args',
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
// PROCESSES: BWA 
// Paired End read alignment 
// Utilizes reference genome
//----------------------------------------

process BWA {
    tag "$sample_id"
    publishDir "${params.outdir}/align/${params.align}",
                pattern: "*.bam", mode: 'copy'
    publishDir "${params.outdir}/align/${params.align}", 
                pattern: "*.sai", mode: 'copy'

    input:
        file ref            
        tuple val(sample_id), path(r1), path(r2)
    
    output:
        file "${sample_id}_1.sai"
        file "${sample_id}_2.sai"
        tuple val(sample_id), path(bam), emit: 'align'
    
    echo true
    
    script:
        r1_sai = "${sample_id}_1.sai"
        r2_sai = "${sample_id}_2.sai"
        align = "${sample_id}_${params.align}_align_pe.sam"
        bam = "${sample_id}_${params.align}_align_pe.bam"
      
    """
    bwa index -a bwtsw $ref;
    bwa aln -t ${task.cpus} -I $ref $r1 > $r1_sai;
    bwa aln -t ${task.cpus} -I $ref $r2 > $r2_sai;
    bwa sampe -P $ref $r2_sai $r2_sai $r1 $r2 > $align \\
    | samtools sort -@${task.cpus} \\
    | samtools view -F4 -b -o $bam \\
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
    publishDir "${params.outdir}/align/${params.align}",
                pattern: "*.bam", mode: 'copy'

    input:
        file ref
        tuple val(sample_id), path(r1), path(r2)

    output:
        file "${sample_id}_${params.align}_align_pe.bam"
        tuple val(sample_id), path(bam), emit: 'align'

    script:
        index_dir = "./index"
        indexes = "./index/index"
        bam = "${sample_id}_bowtie2_align_pe.bam"
    
    """
    sudo mkdir $index_dir;
    sudo chmod 777 $index_dir;
    bowtie2-build $ref $indexes;
    bowtie2 --threads ${task.cpus} \\
            -x $indexes -1 $r1 -2 $r2 \\
    | samtools sort -@${task.cpus} \\
    | samtools view -F4 -b -o $bam 
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
    publishDir "${params.outdir}/align/${params.align}",
                pattern: "*.bam", mode: "copy"
    
    input:
        file ref            
        tuple val(sample_id), path(r1), path(r2)
    
    output:
        file "${sample_id}_${params.align}_align_pe.bam"
        tuple val(sample_id), path(bam), emit: 'align'
    
    script:
        bam = "${sample_id}_${params.align}_align_pe.bam"

    """
    minimap2 -a -t ${task.cpus} \\
        -ax sr $ref $r1 $r2 \\
        | samtools sort -@${task.cpus} \\
        | samtools view -F4 -b -o $bam 
    """
}

//----------------------------------------

//----------------------------------------
// PROCESSES: BAMINDEX
// Index sorted Bam file  
//----------------------------------------

process BAMINDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/align/${params.align}",
                pattern: "*.bai", mode: "copy"

    input: 
        tuple val(sample_id), path(align)

    output:
        file "${sample_id}_${params.align}_align_pe.bai"
        tuple val(sample_id),  path(align), emit: 'align'
    
    script: 
        indexed_bam = "${sample_id}_${params.align}_align_pe.bai"

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
                pattern: "${sample_id}*.vcf", mode: 'copy'
    publishDir "${params.outdir}/variant/freebayes", 
                pattern: "${sample_id}*.txt", mode: 'copy'
    echo true 
    input: 
        file ref
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path(variant), path(ref), 
                emit: 'variant'
    
    script: 
        variant = "${sample_id}.vcf"
        contamination = "${sample_id}_contamination.txt"

    """
    freebayes --gvcf --use-mapping-quality \\
    --genotype-qualities \\
    -f $ref -g 1000 \\
     $bam > $variant 
    """
    //--contamination-estimates $contamination \\ 
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
    bcftools stats -F $ref $variant -v > $stats; 
    
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
        tuple val(sample_id), path(filtered),
                emit: 'variant'
        file "${sample_id}_filtered.vcf"
            
    script: 
        options = "${params.filter_args_str}"
        filtered = "${sample_id}_filtered.vcf"

    """
    bcftools filter $options $variant -o $filtered; 
    
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
    // Read Mapping
    //----------------------------------------
    
    ch_ref = Channel.value(file("${baseDir}/data/ref.fa"))

    if (params.align == 'bowtie2') {
        
        BOWTIE2(ch_ref, FASTP.out.reads)
        BAMINDEX(BOWTIE2.out.align)
    
    } else if (params.align == 'minimap2') {
        
        MINIMAP2(ch_ref, FASTP.out.reads)
        BAMINDEX(MINIMAP2.out.align)
    
    } else {
        
        BWA(ch_ref, FASTP.out.reads)
        BAMINDEX(BWA.out.align)
    }
    
    //----------------------------------------
    // VARIANT CALLING 
    //----------------------------------------
    
    if (params.variant == 'other_variant_calling _tool') {

        FREEBAYES(ch_ref, BAMINDEX.out.align)
    
    } else {
    
        FREEBAYES(ch_ref, BAMINDEX.out.align)
    
    }
    
    //----------------------------------------
    // STATS
    //----------------------------------------

    BCFTOOLS_STATS(FREEBAYES.out.variant)
    
    //----------------------------------------
    // FILTERING 
    //----------------------------------------

    if (params.filter == 'other_filtering_tool') {

        BCFTOOLS_FILTER(FREEBAYES.out.variant)
    
    } else {
        
        BCFTOOLS_FILTER(FREEBAYES.out.variant)
    }
     
    //----------------------------------------
    // Annotate Genomic Variants    
    //----------------------------------------

    if (params.prediction == 'other_annotation_tool') {
        
        SNPEFF(BCFTOOLS_FILTER.out.variant)
    
    } else {

        SNPEFF(BCFTOOLS_FILTER.out.variant)

    }
    //----------------------------------------
    // Consensus building 
    //----------------------------------------
    
    if (params.consensus == 'vcf_consensus') {

        VCF_CONSENSUS(ch_ref, SNPEFF.out.annotation)
    
    } else {

        BCFTOOLS_CONSENSUS(ch_ref, SNPEFF.out.annotation)
    
    }
}

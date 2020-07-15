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
// HELPER FUNCTIONS: PRINT TOOL ARGS
//----------------------------------------
def print_tool_args(tool, tool_args) {
    println(tool_args)
    println("""$tool: ${params[tool]}
    => ${tool_args[0]}: ${params[tool_args[0]]}\n""") 
}
//----------------------------------------

//=============================================================================
// INPUT VALIDATION
//=============================================================================

params.tool_args.each{ k, v ->  print_tool_args(k, v) }

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

//=============================================================================
// PROCESSES
//=============================================================================
include BWA from "./modules/align.nf"
include BOWTIE2 from "./modules/align.nf"
include MINIMAP2 from "./modules/align.nf"

//----------------------------------------
// PROCESSES: BAMINDEX
// Index sorted Bam file  
//----------------------------------------
process BAMINDEX {
    tag "$sample_id"
    publishDir "${align_dir}",
                pattern: "*.bai", mode: "copy"

    echo true 

    input: 
        tuple val(method), val(sample_id), 
              path(align), val(align_dir)

    output:
        file "${sample_id}_${params.align}_align_pe.bai"
        tuple val(method), val(sample_id),  path(align), 
              path(depths), emit: 'align'
    
    script: 
        indexed_bam = "${sample_id}_${params.align}_align_pe.bai"
        depths = "${sample_id}_depths.tsv"

    """
    samtools index $align $indexed_bam;
    samtools depth -a -d 0 $align \\
    | perl -ne 'chomp \$_; print "${sample_id}\t\$_\n"' > $depths;
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
    publishDir "${params.outdir}/variant/${params.variant}/$method",
                pattern: "${sample_id}*.vcf", mode: 'copy'
    publishDir "${params.outdir}/logs/${params.variant}/$method",
                pattern: "${sample_id}.log", mode: "copy"

    input: 
        file ref
        tuple val(method), val(sample_id), 
              path(bam), path(depths)

    output:
        file "${sample_id}.log"
        tuple val(method), val(sample_id), path(variant), path(ref), 
              path(depths), emit: 'variant'
    
    script: 
        variant = "${sample_id}.vcf"
        freebayes_log = "${sample_id}.log"

    """
    freebayes --gvcf --use-mapping-quality \\
    --genotype-qualities \\
    -f $ref -g 1000 $bam > $variant; 
    cat .command.log | tee $freebayes_log;
    """
}
//----------------------------------------

//----------------------------------------
// PROCESSES: BCFTOOLS_STATS
// Collect stats from variant calls 
//----------------------------------------
process BCFTOOLS_STATS {
    tag "$sample_id"
    publishDir "${params.outdir}/variant/stats/$method",
                pattern: "*_stats.txt", mode: 'copy'

    input: 
        tuple val(method), val(sample_id), path(variant), 
              path(ref), path(depths)

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
    publishDir "${params.outdir}/logs/${params.filter}",
                pattern: "${sample_id}.log", mode: "copy"

    echo true 

    input: 
        tuple val(method), val(sample_id), path(variant), 
              path(ref), path(depths)

    output:
        file "${sample_id}_filtered.vcf"
        file "${sample_id}.log"
        tuple val(method), path(filtered), path(ref), path(depths), 
              emit: 'variant'
            
    script: 
        options = "${params.filter_args}"
        filtered = "${sample_id}_filtered.vcf"
        bcftools_filter_log = "${sample_id}.log"
    
    """
    bcftools filter $options $variant -o $filtered; 
    cat .command.log | tee $bcftools_filter_log;
    """
}
//----------------------------------------

//----------------------------------------
// PROCESSES: SNPEFF
// Annotates variants 
// Note: This process is currently deprecated
//----------------------------------------
process SNPEFF {
    tag "$sample_id"
    publishDir "${params.outdir}/variant/prediction",
                pattern: "*_annotation.vcf", mode: 'copy'

    input: 
        tuple val(method), path(variant), path(ref), path(depths)
        tuple val(sample_id), path(r1), path(r2)
    output:
        file "${sample_id}_annotation.vcf"
         tuple val(sample_id), path(variant), path(depths),
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
// PROCESSES: SNIPPY
// Annotates variants 
//----------------------------------------
process SNIPPY {
    tag "$sample_id"
    publishDir "${params.outdir}/variant/snippy/$method",
                pattern: "*snps.*", mode: 'copy'
    publishDir "${params.outdir}/logs/${params.prediction}/$method",
                pattern: "${sample_id}.log", mode: "copy"
    
    input: 
        tuple val(method), path(variant), path(ref), path(depths)
        tuple val(sample_id), path(r1), path(r2)

    output:
        file "${sample_id}.log"
        tuple val(method), val(sample_id), path(variant), path(depths),
                emit: 'annotation' 
    
    script:
        outdir = "${params.outdir}/variant/snippy/$method"
        snippy_log = "${sample_id}.log"
        //travis = (params.travis == 'true')? '--ram 4': ''
    """
    
    snippy --cpus ${task.cpus} --ram 4 \\
    --outdir $outdir \\
    --ref $ref --R1 $r1 --R2 $r2;
    cat .command.log | tee $snippy_log
    
    """
}
//----------------------------------------

//----------------------------------------
// PROCESSES: BCFTOOLS_CONSENSUS
// Build consensus sequence from variants
//----------------------------------------
process BCFTOOLS_CONSENSUS {
    tag "$sample_id"
    publishDir "${params.outdir}/consensus/${params.consensus}/$method",
                pattern: "*.fa", mode: 'copy'
    publishDir "${params.outdir}/logs/${params.consensus}/$method",
                pattern: "${sample_id}.log", mode: "copy"

    input: 
        file ref 
        tuple val(method), val(sample_id), path(variant), path(depths)

    output:
        file "${sample_id}.log"
        file "${sample_id}_consensus.vcf"
       
    script:
        zipped = "${variant}.gz"
        consensus = "${sample_id}_${params.consensus}_consensus.fa"
        bcftools_consensus_log = "${sample_id}.log"
    """
    bgzip $variant;
    tabix -p vcf $zipped;
    bcftools consensus -f $ref $zipped > $consensus;
    cat .command.log | tee $bcftools_consensus_log;
    """
}
//----------------------------------------

//----------------------------------------
// PROCESSES: VCF_CONSENSUS
// Build consensus sequence from variants
// Note: vcf_consensus is the default 
//----------------------------------------
process VCF_CONSENSUS {
    tag "$sample_id"
    publishDir "${params.outdir}/consensus/vcf/$method",
                pattern: "*.fa", mode: 'copy'
    publishDir "${params.outdir}/logs/vcf/$method", 
                pattern: "${sample_id}.log", mode: "copy" 
    input:
        file ref
        tuple val(method), val(sample_id), path(variant), path(depths)

    output:
        file "${sample_id}.log"
        file "${sample_id}_vcf_consensus.fa"

    script:
        consensus = "${sample_id}_vcf_consensus.fa"
        vcf_consensus_log = "${sample_id}.log"

    """
    vcf_consensus_builder \\
        -v $variant \\
        -d $depths \\
        -r $ref \\
        -o $consensus \\
        --sample-name $sample_id;
    cat .command.log | tee $vcf_consensus_log;
    """
}
//----------------------------------------

//=============================================================================
// WORKFLOW DEFINITION 
//=============================================================================
include quality_check from "./modules/quality_check.nf"

//----------------------------------------
// WORKFLOW: variants
// Takes:
//      ref = reference genome
//      align = sorted bam file 
//      
// Main: 
//      1. Indeexes bam file (align)
//      2. Variant calling
//      3. Variant stats
//      4. Filter variants
//
// Emit:
//      variant = Main.4.out.variant 
//
//----------------------------------------
workflow variants {

    take: 
        ref
        align

    main: 
        BAMINDEX(align) 
    
        if (params.variant == 'other_variant_calling _tool') {

            FREEBAYES(ref, BAMINDEX.out.align)
    
        } else {
    
            FREEBAYES(ref, BAMINDEX.out.align)
    
        }
    
        BCFTOOLS_STATS(FREEBAYES.out.variant)
    
        if (params.filter == 'other_filtering_tool') {

            BCFTOOLS_FILTER(FREEBAYES.out.variant)
    
        } else {
        
            BCFTOOLS_FILTER(FREEBAYES.out.variant)
        }
            
    emit:
        variant = BCFTOOLS_FILTER.out.variant

}
//----------------------------------------

//----------------------------------------
// WORKFLOW: consensus
// Takes:
//      variants = variant calls on reads
//      reads = paired end reads
//      ref = reference genome
//      
// Main: 
//      1. Annotation && effect prediction
//      2. Consensus
//
//----------------------------------------
workflow consensus {

    take:
        variants
        reads
        ref

    main: 
       
        if (params.prediction == 'snpeff') {
            //deprecated until further notice        
            SNPEFF(variants, reads)
        } else {

            SNIPPY(variants, reads)
        }
    
        if (params.consensus == 'bcftools') {

            BCFTOOLS_CONSENSUS(ref, SNIPPY.out.annotation)
    
        } else {

            VCF_CONSENSUS(ref, SNIPPY.out.annotation)
    
        }
}
//----------------------------------------

//----------------------------------------
// WORKFLOW: bwa
// Takes:
//      ref = reference genome
//      reads = paired end reads
//      phix = phix genome
//      
// Main: 
//      1. Quality check and read filtering 
//      2. Read mapping using BWA
//      3. Variant calling
//      4. Consensus 
//
//----------------------------------------
workflow bwa {

    take: 
        ref
        reads
        phix

    main:
        quality_check(reads, phix)
        BWA(ref, quality_check.out.reads)
        variants(ref, BWA.out.align)
        consensus(variants.out.variant,
                  reads, ref)
}
//----------------------------------------

//----------------------------------------
// WORKFLOW: bowtie2
// Takes:
//      ref = reference genome
//      reads = paired end reads
//      phix = phix genome
//      
// Main: 
//      1. Quality check and read filtering 
//      2. Read mapping using BOWTIE2
//      3. Variant calling
//      4. Consensus 
//
//----------------------------------------
workflow bowtie2 {

    take: 
        ref
        reads
        phix

    main:
        quality_check(reads, phix)
        BOWTIE2(ref, quality_check.out.reads)
        variants(ref, BOWTIE2.out.align)
        consensus(variants.out.variant, reads, ref)
}
//----------------------------------------

//----------------------------------------
// WORKFLOW: minimap2
// Takes:
//      ref = reference genome
//      reads = paired end reads
//      phix = phix genome
//      
// Main: 
//      1. Quality check and read filtering 
//      2. Read mapping using MINIMAP2
//      3. Variant calling
//      4. Consensus 
//
//----------------------------------------
workflow minimap2 {

    take: 
        ref
        reads
        phix 

    main:
        quality_check(reads, phix)
        MINIMAP2(ref, quality_check.out.reads)
        variants(ref, MINIMAP2.out.align)
        consensus(variants.out.variant,reads, ref)
}

//----------------------------------------

workflow {
    main:    
        ref = Channel.value(file("${params.ref}"))
        
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

        if(params.align == 'bowtie2' || params.align == 'all')
            bowtie2(ref, reads, phix)

        if(params.align == 'minimap2' || params.align == 'all')
            minimap2(ref, reads, phix) 
        
        if(params.align == 'bwa' || params.align == 'all')
            bwa(ref, reads, phix)

}

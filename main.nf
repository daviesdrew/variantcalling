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
// HELPER FUNCTIONS: TOOL ARGS
//----------------------------------------
params.tool_args = [ 'align': 'align_args', 
                     'variant': 'variant_args',
                     'filter': 'filter_args',
                     'prediction': 'prediction_args',
                     'consensus': 'consensus_args' ] 

//----------------------------------------

//=============================================================================
// HELPER FUNCTIONS
//=============================================================================

//----------------------------------------
// HELPER FUNCTIONS: CHECK ARG EXISTENCE
//----------------------------------------
def check_arg_existence(tool, tool_args) {
    params[tool] = (params.containsKey(tool)) ?
                        params[tool] :
                        null
    params[tool_args] = (params.containsKey(tool_args)) ?
                        params[tool_args] :
                        ""
}
//----------------------------------------

//----------------------------------------
// HELPER FUNCTIONS: PRINT TOOL ARGS
//----------------------------------------
def print_tool_args(k, v) {
    println("$k: ${params[k]} && $v: ${params[v]}") 
}
//----------------------------------------

//=============================================================================
// INPUT VALIDATION
//=============================================================================

params.tool_args.each{ k, v -> check_arg_existence(k, v) }
params.tool_args.each{ k, v -> print_tool_args(k, v) }

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

    input:
        file ref            
        tuple val(sample_id), path(r1), path(r2)
    
    output:
        tuple val(sample_id), path(bam), emit: 'align'
    
    script:
        align = "${sample_id}_${params.align}_align_pe.sam"
        bam = "${sample_id}_${params.align}_align_pe.bam"
      
    """
    bwa index -a bwtsw $ref;
    bwa mem -P -t ${task.cpus} $ref $r1 $r2 -o $align;  
    samtools sort $align -@${task.cpus} \\
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
        tuple val(sample_id),  path(align), path(depths), 
              emit: 'align'
    
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
    publishDir "${params.outdir}/variant/freebayes",
                pattern: "${sample_id}*.vcf", mode: 'copy'
    publishDir "${params.outdir}/variant/freebayes", 
                pattern: "${sample_id}*.txt", mode: 'copy'
    echo true 
    input: 
        file ref
        tuple val(sample_id), path(bam), path(depths)

    output:
        tuple val(sample_id), path(variant), path(ref), 
              path(depths), emit: 'variant'
    
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
        tuple val(sample_id), path(variant), 
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

    input: 
        tuple val(sample_id), path(variant), 
              path(ref), path(depths)

    output:
        tuple val(sample_id), path(filtered), path(depths),
                emit: 'variant'
        file "${sample_id}_filtered.vcf"
            
    script: 
        options = "${params.filter_args}"
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
        tuple val(sample_id), path(variant), path(depths)

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
    publishDir "${params.outdir}/variant/prediction",
                pattern: "*_annotation.vcf", mode: 'copy'

    input: 
        tuple val(sample_id), path(variant), path(depths)

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
// PROCESSES: BCFTOOLS_CONSENSUS
// Build consensus sequence from variants
//----------------------------------------
process BCFTOOLS_CONSENSUS {
    tag "$sample_id"
    publishDir "${params.outdir}/consensus/${params.consensus}",
                pattern: "*.fa", mode: 'copy'

    input: 
        file ref 
        tuple val(sample_id), path(variant), path(depths)

    output:
        file "${sample_id}_consensus.vcf"
       
    script:
        zipped = "${variant}.gz"
        consensus = "${sample_id}_${params.consensus}_consensus.fa"

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
                pattern: "*.fa", mode: 'copy'
    echo true 
    
    input:
        file ref
        tuple val(sample_id), path(variant), path(depths)

    output:
        file "${sample_id}_vcf_consensus.fa"

    script:
        consensus = "${sample_id}_vcf_consensus.fa"
        

    """
    vcf_consensus_builder \\
        -v $variant \\
        -d $depths \\
        -r $ref \\
        -o $consensus \\
        --sample-name $sample_id
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
    
    ch_ref = Channel.value(file("${baseDir}/${params.ref}"))

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
    
    if (params.consensus == 'bcftools') {

        BCFTOOLS_CONSENSUS(ch_ref, SNPEFF.out.annotation)
    
    } else {

        VCF_CONSENSUS(ch_ref, SNPEFF.out.annotation)
    
    }
}

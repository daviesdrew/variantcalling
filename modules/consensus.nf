//=============================================================================
// PROCESSES
//=============================================================================

//----------------------------------------
// PROCESSES: SNIPPY
// Annotates and predicts variants
//----------------------------------------
process SNIPPY {
    tag "$sample_id"

    publishDir "${params.outdir}/variant/snippy", 
                pattern: "${file_base}*", mode: "copy"
    publishDir "${params.outdir}/logs/${params.prediction}",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "${logfile}" }
    
    input:
        tuple val(file_base), path(variant), path(ref)
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), val(file_base), 
              path(variant), emit: 'annotation'
        file ".command.log"
    
    script:
        outdir = "${params.outdir}/variant/snippy"
        prefix = "${file_base}"
        logfile = "${file_base}.log"
        travis = (params.travis == true) ? '--ram 4' : ''
        
    """
    snippy --cpus ${task.cpus} $travis \\
    --outdir $outdir --prefix $prefix \\
    --ref $ref --R1 $r1 --R2 $r2;
    """
}
//----------------------------------------

//----------------------------------------
// PROCESSES: BCFTOOLS_CONSENSUS
// Build consensus sequence from variants
//----------------------------------------
process BCFTOOLS_CONSENSUS {
    tag "$sample_id"

    publishDir "${params.outdir}/consensus/bcf",
                pattern: "*.fa", mode: "copy"
    publishDir "${params.outdir}/logs/bcf",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "${logfile}" }
    
    input:
        file ref
        tuple val(sample_id), val(file_base), path(variant)

    output:
        tuple val(sample_id), val(file_base), 
              path(variant), emit: 'consensus'
        file "${consensus}"
        file ".command.log"

    script:
        temp = "temp.vcf"
        zipped = "${variant}.gz"
        consensus = "${file_base}_CONSENSUS.fa"
        logfile = "${file_base}.log"

    """
    cp $variant $temp;
    bgzip $variant;
    tabix -p vcf $zipped;
    bcftools consensus -f $ref $zipped > $consensus;
    cp $temp $variant;
    """
}
//----------------------------------------
// PROCESSES: VCF_CONSENSUS
// Build consensus sequence from variants
//----------------------------------------
process VCF_CONSENSUS {
    tag "$sample_id"

    publishDir "${params.outdir}/consensus/vcf",
                pattern: "*.fa", mode: "copy"
    publishDir "${params.outdir}/logs/vcf",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "${logfile}" }

    input: 
        file ref 
        tuple val(sample_id), val(file_base), path(variant)
        path depths 

    output: 
        tuple val(sample_id), val(file_base), 
              path(variant), emit: 'consensus'
        file "${consensus}"
        file ".command.log"

    script:
        consensus = "${file_base}_CONSENSUS.fa"
        logfile = "{file_base}.log"

    """
    vcf_consensus_builder \\
        -v $variant \\
        -d $depths \\
        -r $ref \\
        -o $consensus \\
        --sample-name $sample_id;
    """
}
//----------------------------------------

//=============================================================================
// WORKFLOW DEFINITION
//=============================================================================

//----------------------------------------
// WORKFLOW: build_consensus
// 
// Takes:
//      ref = reference genome
//      annotation = annotated variant calls
//      depths = depths of reads
//
// Main: 
//      1. Build Consensus
//
// Emit:
//      consensus = consensus sequence
//
//----------------------------------------
workflow build_consensus {
    take:
        ref
        annotation
        depths

    main:
        BCFTOOLS_CONSENSUS(ref, annotation)
        VCF_CONSENSUS(ref, annotation, depths)
    
    emit:
        vcf =  VCF_CONSENSUS.out.consensus 
        bcftools = BCFTOOLS_CONSENSUS.out.consensus
}
//----------------------------------------

//----------------------------------------
// WORKFLOW: consensus
// 
// Takes:
//      variants = variant calls on reads
//      reads = paired end reads
//      ref = reference genome
//
// Main: 
//      1. Annotation && effect prediction 
//      2. Build consensus
//
//----------------------------------------
workflow consensus {

    take:
        variants
        depths
        reads
        ref

    main:

        SNIPPY(variants, reads) 
        build_consensus(ref, SNIPPY.out.annotation, depths)

    emit:
        vcf = build_consensus.out.vcf
        bcf = build_consensus.out.bcftools
}
//----------------------------------------

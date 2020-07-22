//=============================================================================
// PROCESSES
//=============================================================================

//----------------------------------------
// PROCESSES: SNIPPY
// Annotates and predicts variants
//----------------------------------------
process SNIPPY {
    tag "$sample_id"

    publishDir "${params.outdir}/variant/snippy/$method", 
                pattern: "*snps.*", mode: "copy"
    publishDir "${params.outdir}/logs/${params.prediction}/$method",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "snippy_${sample_id}.log" }

    input:
        tuple val(method), path(variant), path(ref), path(depths)
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(method), val(sample_id),
                path(variant), path(depths), emit: 'annotation'
        file ".command.log"
    
    script:
        outdir = "${params.outdir}/variant/snippy/$method"
        
    """
    snippy --cpus ${task.cpus} \\
    --outdir $outdir \\
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

    publishDir "${params.outdir}/consensus/${params.consensus}/$method",
                pattern: "*.fa", mode: "copy"
    publishDir "${params.outdir}/logs/${params.consensus}/$method",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "bcftools_consensus_${sample_id}.log" }
    
    input:
        file ref
        tuple val(method), val(sample_id),
                path(variant), path(depths)

    output:
        file "${sample_id}_bcftools_consensus.fa"
        file ".command.log"

    script:
        zipped = "${variant}.gz"
        consensus = "${sample_id}_bcftools_consensus.fa"

    """
    bgzip $variant;
    tabix -p vcf $zipped;
    bcftools consensus -f $ref $zipped > $consensus;
    """
}
//----------------------------------------
// PROCESSES: VCF_CONSENSUS
// Build consensus sequence from variants
//----------------------------------------
process VCF_CONSENSUS {
    tag "$sample_id"

    publishDir "${params.outdir}/consensus/vcf/$method",
                pattern: "*.fa", mode: "copy"
    publishDir "${params.outdir}/logs/vcf/$method",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "vcf_consensus_${sample_id}.log" }

    input: 
        file ref 
        tuple val(method), val(sample_id),
                path(variant), path(depths) 

    output: 
        file "${sample_id}_vcf_consensus.fa"
        file ".command.log"

    script:
        consensus = "${sample_id}_vcf_consensus.fa"

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
// WORKFLOW: consensus
// 
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
        if (params.prediction == 'different annotation/prediction tool') {

            SNIPPY(variants, reads)

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


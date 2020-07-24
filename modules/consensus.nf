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

    publishDir "${params.outdir}/consensus/${params.consensus}",
                pattern: "*.fa", mode: "copy"
    publishDir "${params.outdir}/logs/${params.consensus}",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "${logfile}" }
    
    input:
        file ref
        tuple val(sample_id), val(file_base), path(variant)

    output:
        file "${consensus}"
        file ".command.log"

    script:
        zipped = "${variant}.gz"
        consensus = "${file_base}_CONSENSUS.fa"
        logfile = "${file_base}.log"

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
        depths
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

            VCF_CONSENSUS(ref, SNIPPY.out.annotation, depths)

        }

}
//----------------------------------------


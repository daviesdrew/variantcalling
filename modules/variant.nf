//=============================================================================
// PROCESSES
//=============================================================================

//----------------------------------------
// PROCESSES: BAMINDEX
// Index sorted bam file
//----------------------------------------
process BAMINDEX {
    tag "$sample_id"

    publishDir "${params.outdir}/align", 
                pattern: "*.bai", mode: "copy"
    
    input:
        tuple val(sample_id), val(file_base), path(align)

    output:
        tuple val(sample_id), val(file_base), path(align), emit: 'align'
        path depths, emit: 'depths'
        file "${indexed_bam}"
        file "${depths}"

    script:
        indexed_bam = "${file_base}.bai"
        depths = "${file_base}_DEPTHS.tsv"

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

    publishDir "${params.outdir}/variant/${params.variant}",
                pattern: "${sample_id}*.vcf", mode: "copy"
    publishDir "${params.outdir}/logs/${params.variant}",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "${logfile}" }
    
    input:
        file ref
        tuple val(sample_id), val(file_base), path(bam)

    output:
        tuple val(sample_id), val(file_base), path(variant), 
              path(ref), emit: 'variant'
        file ".command.log"
    
    script:
        variant = "${file_base}.vcf"
        logfile = "${file_base}.log"

    """
    freebayes --gvcf --use-mapping-quality \\
    --genotype-qualities \\
    -f $ref -g 1000 $bam > $variant;
    """
}
//----------------------------------------

//----------------------------------------
// PROCESSES: BCFTOOLS_STATS
// Collect stats from variant calls
//----------------------------------------
process BCFTOOLS_STATS {
    tag "$sample_id"

    publishDir "${params.outdir}/variant/stats",
                pattern: "*_STATS.txt", mode: "copy"
    input: 
        tuple val(sample_id), val(file_base), path(variant), path(ref) 

    output: 
        file "${stats}"

    script:
        stats = "${file_base}_STATS.txt"

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
                pattern: "*_filtered.txt", mode: "copy"
    publishDir "${params.outdir}/logs/${params.filter}",
                pattern: "${sample_id}.log", 
                mode: "copy",
                saveAs: { file -> "${logfile}" } 

    
    input: 
        tuple val(sample_id), val(file_base), path(variant), path(ref) 

    output: 
        tuple val(file_base), path(filtered), path(ref), emit: 'variant'
        file "${filtered}" 
        file ".command.log"

    script:
        options = "${params.filter_args}"
        filtered = "${file_base}_FILTERED.vcf"
        logfile = "${file_base}.log"

    """
    bcftools filter $options $variant -o $filtered;
    """
}
//----------------------------------------

//=============================================================================
// WORKFLOW DEFINITION
//=============================================================================

//----------------------------------------
// WORKFLOW: variants
// 
// Takes:
//      refs = reference genome
//      align = sorted bam file
//
// Main: 
//      1. Indexes bam file (align)
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
        FREEBAYES(ref, BAMINDEX.out.align)
        BCFTOOLS_STATS(FREEBAYES.out.variant)
        BCFTOOLS_FILTER(FREEBAYES.out.variant)

    emit: 
        depths = BAMINDEX.out.depths
        variant = BCFTOOLS_FILTER.out.variant

}
//----------------------------------------


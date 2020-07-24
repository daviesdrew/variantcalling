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
        tuple val(sample_id), val(method), 
              val(ref_base), path(align)

    output:
        tuple val(sample_id), val(method),
              path(align), path(depths), 
              emit: 'align'
        file "${indexed_bam}"

    script:
        file_base = "{method}_${sample_id}_TO_${ref_base}" 
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

    publishDir "${params.outdir}/variant/${params.variant}/$method",
                pattern: "${sample_id}*.vcf", mode: "copy"
    publishDir "${params.outdir}/logs/${params.variant}/$method",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "freebayes_${sample_id}.log" }
    
    input:
        file ref
        tuple val(sample_id), val(method),
                path(bam), path(depths)

    output:
        tuple val(sample_id), val(method),
                path(variant), path(ref), path(depths), emit: 'variant'
        file ".command.log"

    script:
        variant = "${sample_id}.vcf"

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

    publishDir "${params.outdir}/variant/method",
                pattern: "*_stats.txt", mode: "copy"
    input: 
        tuple val(sample_id), val(method),
                path(variant), path(ref), path(depths) 

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
                pattern: "*_filtered.txt", mode: "copy"
    publishDir "${params.outdir}/logs/${params.filter}",
                pattern: "${sample_id}.log", 
                mode: "copy",
                saveAs: { file -> "bcftools_filter_${sample_id}.log" }

    echo true
    
    input: 
        tuple val(sample_id), val(method),
                path(variant), path(ref), path(depths) 

    output: 
        tuple val(method), path(filtered), path(ref), path(depths),
                emit: 'variant'
        file "${sample_id}_filtered.vcf"
        file ".command.log"

    script:
        options = "${params.filter_args}"
        filtered = "${sample_id}_filtered.vcf"

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

        if (params.variant == 'other variant calling tool') {

            FREEBAYES(ref, BAMINDEX.out.align)

        } else {

            FREEBAYES(ref, BAMINDEX.out.align)

        }

        BCFTOOLS_STATS(FREEBAYES.out.variant)

        if (params.filter == 'other filtering tool') {

            BCFTOOLS_FILTER(FREEBAYES.out.variant)

        } else {

            BCFTOOLS_FILTER(FREEBAYES.out.variant)

        }

    emit: 
        variant = BCFTOOLS_FILTER.out.variant

}
//----------------------------------------


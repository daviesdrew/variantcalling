//=============================================================================
// PROCESSES
//=============================================================================

//----------------------------------------
// PROCESSES: BWA
// Paired End read mapping
// Utilized reference genome
//----------------------------------------
process BWA {
    tag "$sample_id"

    publishDir "${params.outdir}/align/bwa",
                pattern: "*.bam"
    publishDir "${params.outdir}/logs/bwa",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "bwa_mem_${sample_id}.log" } 
                
    input:
        val ref
        tuple sample_id, path(r1), path(r2)

    output:
        tuple val(sample_id), val('bwa'), path(bam), emit: "align"
        file ".command.log"

    script:
        align = "${sample_id}_bwa_align_pe.sam"
        bam = "${sample_id}_bwa_align_pe.bam"

    """
    bwa index -a bwtsw $ref;
    bwa mem -P -t ${task.cpus} $ref $r1 $r2 -o $align;
    samtools sort $align -@${task.cpus} \\
    | samtools view -F4 -b -o $bam;
    """
}
//----------------------------------------

//----------------------------------------
// PROCESSES: BOWTIE2
// Paired End read mapping
// Utilized reference genome
//----------------------------------------
process BOWTIE2 {
    tag "$sample_id"

    publishDir "./index", pattern: "index*bt2*", mode: "copy"
    publishDir "${params.outdir}/align/bowtie2",
                pattern: "*.bam"
    publishDir "${params.outdir}/logs/bowtie2",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "bowtie2_${sample_id}.log" } 
               
    input:
        val ref
        tuple sample_id, path(r1), path(r2)

    output:
        tuple val(sample_id), val('bowtie2'), path(bam), emit: "align"
        file "${sample_id}_bowtie2_align_pe.bam"
        file ".command.log"

    script:
        index_dir = "./index"
        indexes = "./index/index"
        bam = "${sample_id}_bowtie2_align_pe.bam"

    """
    mkdir $index_dir
    chmod 777 $index_dir;
    bowtie2-build $ref $indexes
    bowtie2 --threads ${task.cpus} \\
        -x $indexes -1 $r2 -2 $r2 \\
    | samtools sort -@${task.cpus} \\
    | samtools view -F4 -b -o $bam;
    """
}
//----------------------------------------

//----------------------------------------
// PROCESSES: MINIMAP2
// Paired End read mapping
// Utilized reference genome
//----------------------------------------
process MINIMAP2 {
    tag "$sample_id"

    publishDir "${params.outdir}/align/minimap2",
                pattern: "*.bam"
    publishDir "${params.outdir}/logs/minimap2",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "minimap2_${sample_id}.log" } 

    input:
        val ref
        tuple sample_id, path(r1), path(r2)

    output:
        tuple val(sample_id), val('minimap2'), path(bam), emit: "align"
        file "${sample_id}_minimap2_align_pe.bam"
        file ".command.log"

    script:
        bam = "${sample_id}_minimap2_align_pe.bam"

    """
    minimap2 -a -t ${task.cpus} \\
        -ax sr $ref $r1 $r2 \\
    | samtools sort -@${task.cpus} \\
    | samtools view -F4 -b -o $bam;
    """
}
//----------------------------------------

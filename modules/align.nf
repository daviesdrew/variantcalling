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
                pattern: "*.bam", mode: "copy"
    publishDir "${params.outdir}/logs/bwa",
                pattern: "${sample_id}.log", mode: "copy"
                
    input:
        val ref
        tuple sample_id, path(r1), path(r2)

    output:
        file "${sample_id}.log"
        tuple val('bwa'), val(sample_id), path(bam), val(align_dir), emit: "align"

    script:
        align = "${sample_id}_bwa_align_pe.sam"
        bam = "${sample_id}_bwa_align_pe.bam"
        align_dir = "${params.outdir}/align/bwa"
        bwa_log = "${sample_id}.log"

    """
    bwa index -a bwtsw $ref;
    bwa mem -P -t ${task.cpus} $ref $r1 $r2 -o $align;
    samtools sort $align -@${task.cpus} \\
    | samtools view -F4 -b -o $bam;
    cat .command.log | tee $bwa_log;
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
                pattern: "*.bam", mode: "copy"
    publishDir "${params.outdir}/logs/bowtie2",
                pattern: "${sample_id}.log", mode: "copy"
                
    input:
        val ref
        tuple sample_id, path(r1), path(r2)

    output:
        file "${sample_id}.log"
        file "${sample_id}_bowtie2_align_pe.bam"
        tuple val('bowtie2'), val(sample_id), path(bam), 
                val(align_dir), emit: "align"

    script:
        index_dir = "./index"
        indexes = "./index/index"
        bam = "${sample_id}_bowtie2_align_pe.bam"
        align_dir = "${params.outdir}/align/bowtie2"
        bowtie2_log = "${sample_id}.log"

    """
    mkdir $index_dir
    chmod 777 $index_dir;
    bowtie2-build $ref $indexes
    bowtie2 --threads ${task.cpus} \\
        -x $indexes -1 $r2 -2 $r2 \\
    | samtools sort -@${task.cpus} \\
    | samtools view -F4 -b -o $bam;
    cat .command.log | tee $bowtie2_log;
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
                pattern: "*.bam", mode: "copy"
    publishDir "${params.outdir}/logs/minimap2",
                pattern: "${sample_id}.log", mode: "copy"
                
    input:
        val ref
        tuple sample_id, path(r1), path(r2)

    output:
        file "${sample_id}_minimap2_align_pe.bam"
        file "${sample_id}.log"
        tuple val('minimap2'), val(sample_id), path(bam), val(align_dir), emit: "align"

    script:
        bam = "${sample_id}_minimap2_align_pe.bam"
        align_dir = "${params.outdir}/align/minimap2"
        minimap2_log = "${sample_id}.log"

    """
    minimap2 -a -t ${task.cpus} \\
        -ax sr $ref $r1 $r2 \\
    | samtools sort -@${task.cpus} \\
    | samtools view -F4 -b -o $bam;
    cat .command.log | tee $minimap2_log;
    """
}
//----------------------------------------

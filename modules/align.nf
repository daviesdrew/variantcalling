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

    publishDir "${params.outdir}/align",
                pattern: "*.bam"
    publishDir "${params.outdir}/logs",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "${logfile}" } 
    
    input:
        file ref
        tuple sample_id, path(r1), path(r2)

    output:
        tuple val(sample_id), val(file_base), path(bam), emit: "align"
        file ".command.log"

    script:
        method = "BWA"
        ref_base = "${ref.getBaseName()}"
        file_base = "BWA_${sample_id}_TO_${ref_base}" 
        bam = "${file_base}.bam"
        logfile = "${file_base}.log" 

    """
    bwa index -a bwtsw $ref;
    bwa mem -t ${task.cpus} $ref $r1 $r2 \\
    | samtools sort -@${task.cpus} \\
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

    publishDir "./index", pattern: "index*bt2*", 
                mode: "copy"
    publishDir "${params.outdir}/align",
                pattern: "*.bam"
    publishDir "${params.outdir}/logs",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "${logfile}" }
               
    input:
        file ref
        tuple sample_id, path(r1), path(r2)

    output:
        tuple val(sample_id), val(file_base), path(bam), emit: "align"
        file ".command.log"

    script:
        index_dir = "./index"
        indexes = "./index/index"
        method = "BOWTIE2"
        ref_base = "${ref.getBaseName()}"
        file_base = "BOWTIE2_${sample_id}_TO_${ref_base}"
        bam = "${file_base}.bam"
        logfile = "${file_base}.log"

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

    publishDir "${params.outdir}/align",
                pattern: "*.bam"
    publishDir "${params.outdir}/logs",
                pattern: ".command.log", 
                mode: "copy",
                saveAs: { file -> "${logfile}" }

    input:
        file ref
        tuple sample_id, path(r1), path(r2)

    output:
        tuple val(sample_id), val(file_base), path(bam), emit: "align"
        file ".command.log"

    script:
        method = 'MINIMAP2'
        ref_base = "${ref.getBaseName()}"
        file_base = "${method}_${sample_id}_TO_${ref_base}"
        bam = "${file_base}.bam"
        logfile = "${file_base}.log"
    """
    minimap2 -a -t ${task.cpus} \\
        -ax sr $ref $r1 $r2 \\
    | samtools sort -@${task.cpus} \\
    | samtools view -F4 -b -o $bam;
    """
}
//----------------------------------------

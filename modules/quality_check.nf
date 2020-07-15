//=============================================================================
// PROCESSES
//=============================================================================

//----------------------------------------
// PROCESSES: REMOVE_PHIX
// Remove any Coliphage phi-X714 reads using bbduk
// bbduk and options explained above script call
//----------------------------------------
process REMOVE_PHIX {
    tag "$sample_id"

    publishDir "${params.outdir}/qc/phix_removal",
                pattern: "*.txt", mode: "copy"
    publishDir "${params.outdir}/reads/phix_removed",
                pattern: "*.fastq.qz"
                
    input:
        file phix
        tuple sample_id, path(r1), path(r2)

    output:
        tuple sample_id, path(r1_out), path(r2_out), emit: "reads"
        tuple sample_id, path(stats), emit: "stats"

    script:
        r1_out = "${sample_id}_1.phix_removed.fastq.qz"
        r2_out = "${sample_id}_2.phix_removed.fastq.qz"
        stats = "${sample_id}-remove_phix-stats.txt"

    """
    bbduk.sh \\
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
// Adapter trimming and quality filtering using fastp
// By default, only up to 40% of bases can be below the min Phred score base
// quality threshold of Q15
// fastp and options explained above script call
//----------------------------------------
process FASTP {
    tag "$sample_id"

    publishDir "${params.outdir}/fastp/html",
                pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/fastp/json",
                pattern: "*.json", mode: "copy"
    publishDir "${params.outdir}/reads/fastp",
                pattern: "*.fastp.fastq.qz"
                
    input:
        tuple sample_id, path(r1), path(r2)

    output:
        tuple sample_id, path(r1_out), path(r2_out), emit: "reads"
        tuple sample_id, path(html_report), path(json_report), emit: "report"

    script:
        r1_out = "${sample_id}_1.fastp.fastq.qz"
        r2_out = "${sample_id}_2.fastp.fastq.qz"
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
            mode: "copy",
            saveAs: { filename -> 
                        filename.indexOf(".zip") > 0 
                            ? "zips/$filename"
                            : "$filename"
            }
    
    input: 
        tuple val(sample_id), path(r1), path(r2)

    output: 
        file "*_fastqc.{zip,html}"

    """
    fastqc -q $r1 $r2
    """
}
//----------------------------------------

//=============================================================================
// WORKFLOW DEFINITION
//=============================================================================

//----------------------------------------
// WORKFLOW: quality_check
// 
// Takes:
//      reads = paired end reads
//      phix = phix genome
//
// Main: 
//      1. Filter using phix
//      2. Trim paired end reads
//      3. Check quality trimmed pe reads
//
// Emit:
//      reads = Trimmed paired end reads
//
//----------------------------------------
workflow quality_check {
    
    take: 
        reads
        phix

    main: 
        REMOVE_PHIX(phix, reads)
        FASTP(REMOVE_PHIX.out.reads)
        fastqc_reports = FASTQC(FASTP.out.reads)

    emit:
        reads = FASTP.out.reads
}
//----------------------------------------


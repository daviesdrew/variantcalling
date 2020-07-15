//=============================================================================
// WORKFLOW DEFINITION
//=============================================================================
include BWA from "./align.nf"
include BOWTIE2 from "./align.nf"
include MINIMAP2 from "./align.nf"

include quality_check from "./quality_check.nf"
include variants from "./variant.nf"
include consensus from "./consensus.nf"

//----------------------------------------
// WORKFLOW: bwa
// 
// Takes:
//      ref = reference genome
//      reads = paired end reads
//      phix = phix genome
//
// Main: 
//      1. Quality check and read filtering
//      2. Read mapping using BWA
//      3. Variant calling
//      4. Consensus
//
//----------------------------------------
workflow bwa {

    take:
        ref
        reads
        phix 

    main:
        quality_check(reads, phix)
        BWA(ref, quality_check.out.reads)
        variants(ref, BWA.out.align)
        consensus(variants.out.variant,
                  reads, ref)
}
//----------------------------------------

//----------------------------------------
// WORKFLOW: bowtie2
// 
// Takes:
//      ref = reference genome
//      reads = paired end reads
//      phix = phix genome
//
// Main: 
//      1. Quality check and read filtering
//      2. Read mapping using BOWTIE2
//      3. Variant calling
//      4. Consensus
//
//----------------------------------------
workflow bowtie2 {

    take:
        ref
        reads
        phix 

    main:
        quality_check(reads, phix)
        BOWTIE2(ref, quality_check.out.reads)
        variants(ref, BOWTIE2.out.align)
        consensus(variants.out.variant,
                  reads, ref)
}
//----------------------------------------

//----------------------------------------
// WORKFLOW: minimap2
// 
// Takes:
//      ref = reference genome
//      reads = paired end reads
//      phix = phix genome
//
// Main: 
//      1. Quality check and read filtering
//      2. Read mapping using MINIMAP2
//      3. Variant calling
//      4. Consensus
//
//----------------------------------------
workflow minimap2 {

    take:
        ref
        reads
        phix 

    main:
        quality_check(reads, phix)
        MINIMAP2(ref, quality_check.out.reads)
        variants(ref, MINIMAP2.out.align)
        consensus(variants.out.variant,
                  reads, ref)
}
//----------------------------------------


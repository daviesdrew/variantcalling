//=============================================================================
// WORKFLOW DEFINITION
//=============================================================================
include { BWA } from "./align.nf"
include { BOWTIE2 } from "./align.nf"
include { MINIMAP2 } from "./align.nf"
include { BOKEH_PLOT } from "./plot.nf"

include { variants } from "./variant.nf"
include { consensus } from "./consensus.nf"

//----------------------------------------
// WORKFLOW: bwa
// 
// Takes:
//      ref = reference genome
//      reads = paired end reads
//
// Main: 
//      1. Read mapping using BWA
//      2. Variant calling
//      3. Consensus
//
//----------------------------------------
workflow bwa {

    take:
        ref
        reads

    main:
        BWA(ref, reads)
        variants(ref, BWA.out.align)
        consensus(variants.out.variant,
                  variants.out.depths, 
                  reads, ref)    
        BOKEH_PLOT(consensus.out.consensus)
}        
//----------------------------------------

//----------------------------------------
// WORKFLOW: bowtie2
// 
// Takes:
//      ref = reference genome
//      reads = paired end reads
//
// Main: 
//      1. Read mapping using BOWTIE2
//      2. Variant calling
//      3. Consensus
//
//----------------------------------------
workflow bowtie2 {

    take:
        ref
        reads

    main:
        BOWTIE2(ref, reads)
        variants(ref, BOWTIE2.out.align)
        consensus(variants.out.variant,
                  variants.out.depths, 
                  reads, ref)
        BOKEH_PLOT(consensus.out.consensus)
}
//----------------------------------------

//----------------------------------------
// WORKFLOW: minimap2
// 
// Takes:
//      ref = reference genome
//      reads = paired end reads
//
// Main: 
//      1. Read mapping using MINIMAP2
//      2. Variant calling
//      3. Consensus
//
//----------------------------------------
workflow minimap2 {

    take:
        ref
        reads

    main:
        MINIMAP2(ref, reads)
        variants(ref, MINIMAP2.out.align)
        consensus(variants.out.variant,
                  variants.out.depths,
                  reads, ref)
        BOKEH_PLOT(consensus.out.consensus)
}
//----------------------------------------


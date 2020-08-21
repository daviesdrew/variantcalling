//=============================================================================
// PROCESSES 
//=============================================================================

//----------------------------------------
// PROCESSES: BOKEHPLOT
// Plot pipeline data
//----------------------------------------
process BOKEHPLOT {
    tag "$sample_id"

    input:
        tuple val(sample_id), val(file_base), path(variant)

    echo true

    """ 
    echo "$sample_id"
    echo "$file_base"
    echo "$variant"
    """
}
//----------------------------------------

//=============================================================================
// WORKFLOW DEFINITION 
//=============================================================================

//----------------------------------------
// WORKFLOW: vcf_plot
//----------------------------------------
workflow vcf_plot {

    take:
        data

    main:
        plot(data)
}
//----------------------------------------

//----------------------------------------
// WORKFLOW: bcf_plot
//----------------------------------------
workflow bcf_plot {

    take:
        data

    main:
        plot(data)
}
//----------------------------------------

//----------------------------------------
// WORKFLOW: plot
// 
// Takes: 
//      x
//      y
//      z
//
// Main:
//      1.
//      2. 
//      3.
//      4.
//
//----------------------------------------
workflow plot {

    take:
        data

    main:
        BOKEHPLOT(data)
}
//----------------------------------------



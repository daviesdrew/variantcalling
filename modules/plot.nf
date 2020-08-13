process BOKEH_PLOT {
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


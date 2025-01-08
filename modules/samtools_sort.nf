nextflow.enable.dsl = 2

process SAMTOOLS_SORT_BYNAME {
    tag "${samfile.getSimpleName()}"
    label 'big_mem'
    publishDir 'results', mode: 'copy', overwrite: true
    
    input: 
    tuple val(sample_id), path(samfile)
    
    output:
    tuple val(sample_id), path("${samfile.getSimpleName()}_byname.sam")
    
    script:
    
    """
    samtools sort -n ${samfile} -O sam -o ${samfile.getSimpleName()}_byname.sam
    """
}

process SAMTOOLS_SORT_BYREF {
    tag "${samfile.getSimpleName()}"
    label 'big_mem'
    publishDir 'results', mode: 'copy', overwrite: true
    maxForks 10
    
    input: 
    tuple val(sample_id), path(samfile)
    
    output:
    tuple val(sample_id), path("${sample_id}_byref.sam")
    
    script:
    
    """
    samtools sort ${samfile} -O sam -o ${sample_id}_byref.sam
    """
}

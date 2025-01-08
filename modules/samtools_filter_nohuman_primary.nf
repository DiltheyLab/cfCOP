nextflow.enable.dsl = 2

process SAMTOOLS_FILTER_NOHUM {
    tag "${sample_id}"
    label 'big_mem'
        
    input: 
    tuple val(sample_id), path(samfile)
    
    output:
    tuple val(sample_id), path("${sample_id}_nohum_primary.sam")
    
    script:
    
    """
    awk '\$3 !~ "sp9606as"' ${samfile} > ${sample_id}_nohum.sam
    samtools view -h -F 2308 ${sample_id}_nohum.sam -o ${sample_id}_nohum_primary.sam
    """
}
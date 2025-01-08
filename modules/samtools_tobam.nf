nextflow.enable.dsl = 2

process SAMTOOLS_TOBAM {
    tag "${samfile.getSimpleName()}"
    label 'big_mem' // not sure about that
    
    input: 
    tuple val(sample), path(samfile)
    
    output:
    tuple val(sample), path("${samfile.getSimpleName()}.bam")
    
    script:
    
    """
    samtools view -@ ${params.threads_samtools} -b ${samfile} -o ${samfile.getSimpleName()}.bam
    """
}

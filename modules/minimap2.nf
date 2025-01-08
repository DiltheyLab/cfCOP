nextflow.enable.dsl = 2

process MINIMAP2_NANO {
    tag "${sample.id}"
    label 'big_mem'
    memory '300 GB'
    maxForks 1
    //publishDir 'results', mode: 'copy', overwrite: true
    //storeDir 'results'
    
    input: 
    val sample
    
    output:
    tuple val(sample), path("${sample.id}_mm2.sam")
    
    when:
        sample.platform == "nanopore" &&
        sample.mapper == "minimap2" &&
        ! sample.sorted_sam &&
        ! sample.sam
    
    script:
    """
    minimap2 -x map-ont -t ${params.threads} -a -o ${sample.id}_mm2.sam ${params.minimap2_index} ${sample.fastq1}
    """

}

process MINIMAP2_ILLU {
    tag "${sample.id}"
    label 'big_mem'
    memory '300 GB'
    maxForks 1
    //publishDir 'results', mode: 'copy', overwrite: true
    //storeDir 'results'
    
    input: 
    val sample
    
    output:
    tuple val(sample), path("${sample.id}_mm2.sam")
    
    when:
        sample.platform == "illumina" &&
        sample.mapper == "minimap2" &&
        ! sample.sorted_sam &&
        ! sample.sam
    
    script:
    
    if (sample.fastq2){ 
    
    """
    minimap2 -x sr -t ${params.threads} -a ${params.minimap2_index} ${sample.fastq1} \\
    | samtools view -bS -@ $task.cpus -o ${sample.id}_fwd.bam -
    minimap2 -x sr -t ${params.threads} -a ${params.minimap2_index} ${sample.fastq2} \\
    | samtools view -bS -@ $task.cpus -o ${sample.id}_rev.bam -
    samtools merge -n ${sample.id}_mm2.bam ${sample.id}_fwd.bam ${sample.id}_rev.bam
    samtools view -h ${sample.id}_mm2.bam > ${sample.id}_mm2.sam
    rm ${sample.id}_*.bam
    """
    
    } else {
    
    """
    minimap2 -x sr -t ${params.threads} -a -o ${sample.id}_mm2.sam ${params.minimap2_index} ${sample.fastq1}
    """
    
    }

}


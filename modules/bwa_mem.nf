nextflow.enable.dsl = 2

process BWA_MEM_NANO {
    tag "${sample.id}"
    label 'big_mem'
	maxForks 1
    //publishDir 'results', mode: 'copy', overwrite: true
    //storeDir 'results'
    
    input: 
    val sample
	
    output:
    tuple val(sample), file("${sample.id}_abwa.sam")
    
    when:
        sample.platform == "nanopore" &&
        sample.mapper == "bwa" &&
        ! sample.sorted_sam &&
        ! sample.sam
    
    script:
    """
    bwa mem -t ${params.threads} -x ont2d -a ${params.bwa_index} ${sample.fastq1} > ${sample.id}_abwa.sam
    """

}


process BWA_MEM_ILLU {
    tag "${sample.id}"
    label 'big_mem'
	maxForks 1
    //publishDir 'results', mode: 'copy', overwrite: true
    //storeDir 'results'
    
    input: 
    val sample
    
    output:
    tuple val(sample), path("${sample.id}_abwa.sam")
    
    when:
        sample.platform == "illumina" &&
        sample.mapper == "bwa" &&
        ! sample.sorted_sam &&
        ! sample.samless
    
    script:
    
    if (sample.fastq2){ 
    
    """
    bwa mem -t ${params.threads} -a ${params.bwa_index} ${sample.fastq1} \\
    | samtools view -bS -@ $task.cpus -o ${sample.id}_fwd.bam -
    bwa mem -t ${params.threads} -a ${params.bwa_index} ${sample.fastq2} \\
    | samtools view -bS -@ $task.cpus -o ${sample.id}_rev.bam -
    samtools merge -n ${sample.id}_abwa.bam ${sample.id}_fwd.bam ${sample.id}_rev.bam
    samtools view -h ${sample.id}_abwa.bam > ${sample.id}_abwa.sam
    rm ${sample.id}_*.bam
    """
    
    } else {
    
    """
    bwa mem -t ${params.threads} -a ${params.bwa_index} ${sample.fastq1} > ${sample.id}_abwa.sam
    """
    
    }

}


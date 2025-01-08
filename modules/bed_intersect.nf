nextflow.enable.dsl = 2

process BEDTOOLS_INTERSECT {
    tag "${filter.getSimpleName()}_${bamfile.getSimpleName()}"
    //publishDir 'results', mode: 'copy', overwrite: true // TEMPORARY! only to test!!! later remove
    errorStrategy = { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries = 3
    memory = { 10.GB * task.attempt }
    //maxForks 5
    
    input: 
    tuple val(sample), path(bamfile), val(filtername), path(filter)
    
    output:
    tuple val(sample), val(filtername), path("${sample.map_id}_${filtername}.txt")
    
    script:
    
    """
    bedtools intersect -wo -bed -a ${bamfile} -b ${filter} > ${sample.map_id}_${filtername}.bed
    less ${sample.map_id}_${filtername}.bed | cut -f 4 | sort | uniq > ${sample.map_id}_${filtername}.txt
    """
}

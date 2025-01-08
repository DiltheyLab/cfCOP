nextflow.enable.dsl = 2

process MAKE_READSTATS {
    tag "${sample.map_id}"
    publishDir 'results', mode: 'copy', overwrite: true // TEMPORARY! only to test!!! later remove
    errorStrategy = { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries = 3
    memory = { 3.GB * task.attempt }
    
    input: 
    tuple val(sample), path(samfile), path(compl_filter), path(hom_filter), path(plasmid_filter)
    path taxonomy
    
    output:
    tuple val(sample), path("${sample.map_id}.readstat")
    
    script:
    
    """
    readstats_from_sam.py ${samfile} ${sample.map_id}.readstat ${taxonomy} --complexity ${compl_filter} --homology ${hom_filter} --plasmid ${plasmid_filter}
    """
}

nextflow.enable.dsl = 2

process MAKE_TAXSTATS {
    tag "${sample.map_id}"
    publishDir 'results', mode: 'copy', overwrite: true // TEMPORARY! only to test!!! later remove
    errorStrategy = { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries = 1
    memory = { 50.GB * task.attempt }
    
    input: 
    tuple val(sample), path(readstat), path(metagen_sam), val(taxrank), val(k), path(taxonomy), path(db_fasta)
    
    output:
    tuple val(sample), val(taxrank), path("${sample.map_id}_${taxrank}.taxstat")
    
    script:
    
    """
    from_readstats_to_taxstats_kmers.py ${readstat} ${taxrank} ${metagen_sam} ${k} ${sample.map_id}_${taxrank}.taxstat ${taxonomy} ${db_fasta}
    """

}
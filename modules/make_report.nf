nextflow.enable.dsl = 2

process MAKE_REPORT {
    tag "${sample.map_id}"
    publishDir 'results', mode: 'copy', overwrite: true // TEMPORARY! only to test!!! later remove
    errorStrategy = { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries = 1
    memory = { 50.GB * task.attempt }
    
    input: 
    tuple val(sample), path(taxstat_spec), path(taxstat_gen), path(taxstat_fam), val(kmer_cutoff), val(count_cutoff), val(fraction_cutoff)
    
    output:
    tuple val(sample), path("${sample.map_id}.report")
    
    script:
    
    """
    family_genus_species_report.py ${taxstat_spec} ${taxstat_gen} ${taxstat_fam} ${sample.map_id}.report  --kmer_cutoff ${kmer_cutoff} --count_cutoff ${count_cutoff} --fraction_cutoff ${fraction_cutoff}
    """

}
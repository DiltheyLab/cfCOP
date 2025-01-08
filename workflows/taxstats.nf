nextflow.enable.dsl = 2

include { SAMTOOLS_SORT_BYREF } from '../modules/samtools_sort'
include { SAMTOOLS_FILTER_NOHUM } from '../modules/samtools_filter_nohuman_primary'
include { MAKE_TAXSTATS } from '../modules/make_taxstats'


workflow TAXSTATS {
    
    take:
        ch_sam
        ch_readstat
        ch_ranks
    
    main:
        ch_sam.map{ sample, sam -> [ sample.map_id, sam ]}.set {named_sam}
        ch_sam.map{ sample, sam -> [ sample.map_id, sample ]}.set {named_sample}
        ch_readstat.map{ sample, readstat -> [ sample.map_id, readstat ]}.set {named_readstat}
        
        // creating metagenomic sams
        ch_nohum_sam = SAMTOOLS_FILTER_NOHUM(named_sam)
        ch_nohum_sorted_sam = SAMTOOLS_SORT_BYREF(ch_nohum_sam)
        
        
        ch_readstat_sam_grouped = named_sample.concat(named_readstat, ch_nohum_sorted_sam).groupTuple()
        
        ch_readstat_sam_grouped.map{ sample, remainder -> remainder}.set {ch_readstat_sam}
        
        taxstat_input = ch_readstat_sam.combine(ch_ranks).combine(Channel.of([params.kmer_size, params.taxonomy, params.database_fasta]))
        
        ch_taxstat_out = MAKE_TAXSTATS(taxstat_input) 
        
            
    emit:
        ch_taxstat_out
        
}


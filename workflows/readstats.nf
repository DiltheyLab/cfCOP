nextflow.enable.dsl = 2

include { SAMTOOLS_TOBAM } from '../modules/samtools_tobam'
include { BEDTOOLS_INTERSECT } from '../modules/bed_intersect'
include { MAKE_READSTATS } from '../modules/make_readstats'

process DUMMY {
    input: 
    tuple val(sample), path(bamfile), val(filtername), path(filter)
    
    output:
    tuple val(sample), val(filtername), path("${sample.map_id}_${filtername}.txt")
    
    script:
    
    """
    echo bla
    """
}



workflow READSTATS {
    
    take:
        ch_sam
        ch_filter
    
    main:
    
        ch_bam = SAMTOOLS_TOBAM(ch_sam)
        ch_bam_filter = ch_bam.combine(ch_filter)
        ch_filtered = BEDTOOLS_INTERSECT(ch_bam_filter)
        
        ch_filtered.map{ sample, filtername, filtered -> [filtername, [sample.map_id, filtered] ]}.set {named_filtered}
    
        named_filtered.branch {
                complexity: it[0] == "complexity"
                homology: it[0] == "homology"
                plasmid: it[0] == "plasmid"
                }
                .set {branched}
        
        complexity = branched.complexity.map {a,b -> b}
        homology = branched.homology.map {a,b -> b}
        plasmid = branched.plasmid.map {a,b -> b}
        
        ch_sam.map{ sample, sam -> [ sample.map_id, sam ]}.set {named_sam}
        ch_sam.map{ sample, sam -> [ sample.map_id, sample ]}.set {named_sample}
        
        // for now all three filters must be present
        // later on alternative processes can be written for only one or two or no filters
        grouped_readstat_input = named_sample.concat(named_sam, complexity, homology, plasmid).groupTuple()
        grouped_readstat_input.map{a,b -> b}.set{ch_readstat_in}
        
        ch_readstat_out = MAKE_READSTATS(ch_readstat_in, params.taxonomy)
        
            
    emit:
        ch_readstat_out
        
}

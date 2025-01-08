nextflow.enable.dsl = 2

include { MAKE_REPORT } from '../modules/make_report'

workflow REPORT {
    
    take:
        ch_taxstats
    
    main:
        ch_taxstats
		.groupTuple()
		.map { sample, taxranks, paths ->
        // Define the desired order of tax ranks
        def desiredOrder = ["species", "genus", "family"]

        // Create a map of taxrank to path
        def taxrankToPath = [taxranks, paths].transpose().collectEntries()

        // Reorder paths based on the desired tax rank order
        def orderedPaths = desiredOrder.collect { taxrankToPath[it] }

        // Return the updated entry with ordered paths
        [sample, orderedPaths]
		}
		.map { sample, nestedList ->
        [sample] + nestedList // Combine sample with the flattened list
		}
		.set { ch_taxstats_in }
		
		ch_report_in = ch_taxstats_in.combine(Channel.of([params.kmer_cutoff, params.count_cutoff, params.fraction_cutoff]))
		
		ch_report_out = MAKE_REPORT(ch_report_in)
    
    emit:
        ch_report_out
        
}

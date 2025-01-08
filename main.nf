nextflow.enable.dsl = 2

include { CHECK_INPUT } from './modules/check_input'
include { path_field_valid } from './modules/check_input'
include { CHECK_DATABASES } from './modules/check_databases'

include { MAP } from './workflows/mapping'

include { READSTATS } from './workflows/readstats'
include { TAXSTATS } from './workflows/taxstats'
include { REPORT } from './workflows/report'

workflow {
	/*
	if (params.sample_sheet) {
        samples = file(params.sample_sheet)
        println "Using sample sheet: $samples"
    } else if (params.input_file) {
        samples = file(params.input_file)
        println "Using input files: $samples"
    } else {
        error "You must provide either a sample sheet (--samplesheet) or input files (--input)"
    }
	*/
	
    samples = CHECK_INPUT(params.input_files, params.samplesheet)
	//samples.view()
	//CHECK_DATABASES(samples)
	
    // ch_sam consists of tuples with sample dictionary and path to sorted samfile
    ch_sam = MAP(samples)
    
    // here is how to use map to transform channels
    ch_sam.map { sample, sam -> [ sample.map_id, sam ]}.groupTuple()
    
    //choosing filters
    ch_filters = Channel.empty()
    
    if (path_field_valid(params.complexity_filter)){
        compl = channel.of(["complexity", file(params.complexity_filter)])
        ch_filters = ch_filters.concat(compl)
        }
    if (path_field_valid(params.homology_filter)){
        hom = channel.of(["homology", file(params.homology_filter)])
        ch_filters = ch_filters.concat(hom)
        }
    if (path_field_valid(params.plasmid_filter)){
        plasm = channel.of(["plasmid", file(params.plasmid_filter)])
        ch_filters = ch_filters.concat(plasm)
        }
    
    ch_readstats = READSTATS(ch_sam, ch_filters)
    
    ch_ranks = Channel.from(params.taxrank)
    
    ch_taxstats = TAXSTATS(ch_sam, ch_readstats, ch_ranks)
	ch_report = REPORT(ch_taxstats)
	//ch_report.view()
	}

nextflow.enable.dsl = 2

// As an inspiration I used:
// https://github.com/nf-core/viralrecon/blob/master/bin/check_samplesheet.py
// https://github.com/nf-core/viralrecon/blob/master/modules/local/samplesheet_check.nf
// https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/input_check.nf


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Sample sheet format: 
    
    1. Unique sample ID (e.g. Project_sampleID_material/comment etc...)
    2. Sequencing platform (illumina/nanopore)
    3. mapper
    4. Fastq 1
    5. Fastq 2 (can be empty if Illumina single end or nanopore)
    6. ggf samfile
    
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def path_field_valid(String somepath){
    if (!somepath) {
        return false
    }
    if (file(somepath).exists()) {
        return true
    }
    else { 
        return false }

}


def create_fastq_channels(LinkedHashMap row){
        // create an empty map
        def mapperlist = ["bwa", "minimap2"]
        def sample = [:]
        
        if (!row.id){
            exit 1, "ERROR: Please check input samplesheet -> Missing sample unique ID\n${row.id}"
        }
        if (!row.platform){
            exit 1, "ERROR: Please check input samplesheet -> Missing sequencing platform\n${row.platform}"
        }
        
        sample.id = row.id
        sample.platform = row.platform
        
        // test if fastq1 exists
        if (!path_field_valid(row.fastq1)) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq1}"
        }
        sample.fastq1 = file(row.fastq1)
        
        // test if fastq2 exists
        if (path_field_valid(row.fastq2)) {
            sample.fastq2 = file(row.fastq2)
        }
        
        // if no or invalid mapper specified, the sample is mapped with default params.mapper
        if (mapperlist.contains(row.mapper)) {
            sample.mapper = row.mapper
            } else {
                sample.mapper = params.mapper
            }
            
        sample.map_id = sample.id + "_" + sample.mapper  
        
        if (path_field_valid(row.sorted_sam)) {
            sample.sorted_sam = row.sorted_sam
            }
        
        // return map with id, platform, fastq path, ggf mapper
        return sample
        
}

def create_multi_file_channel(file_name) {
    // Split the input by comma and trim whitespace
    def files = file_name.toString().split(",").collect { it.trim() }
    
    // Map each file to a LinkedHashMap structure
    return files.collect { infile ->
        def clean_name = infile.toString().replaceAll(/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$/, "")
                                        .replaceAll(/[^a-zA-Z0-9]/, "_")

        // Create LinkedHashMap for each file
        return new LinkedHashMap([
            id     : clean_name,            // Sample ID
            platform: params.platform,
			mapper: params.mapper,
			fastq1 : infile,                // FastQ file (single-end or first pair)
			fastq2 : null                   // Placeholder for paired-end files
        ])
    }
}

workflow CHECK_INPUT {
    take:
    samplesheet  
    input_files  

    main:
    samples = params.samplesheet ? 
        Channel.fromPath(params.samplesheet)
            .splitCsv(header:true, sep:"\t")
            .map { create_fastq_channels(it) } :
        params.input_files ?
        Channel.from(params.input_files)
            .flatMap { create_multi_file_channel(it) } // Creates LinkedHashMap objects
            .map { create_fastq_channels(it) } :       // Pass each LinkedHashMap to the function

        error("You must provide either a sample sheet or input files.")
    
    emit:
    samples 
}

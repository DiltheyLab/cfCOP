nextflow.enable.dsl = 2
/*
def verify_bwa_database () {
	if (!params.bwa_index) {
			throw new IllegalArgumentException("Please specify path to BWA database. Parameter bwa_index is not specified or is null.")
		}
		
	//println "Workflow parameters: ${params.bwa_index}"
	
	ch_bwa_db = Channel.fromPath("${params.bwa_index}.{amb,ann,bwt,pac,sa}").toList()
	//ch_bwa_db.view()
	
	ch_bwa_db.then { fileList ->
		if (fileList == null || fileList.isEmpty()) {
			throw new IllegalArgumentException("Workflow terminated! No valid files found for the specified BWA database path: ${params.bwa_index}")
		}
		if (fileList.size() != 5) {
			throw new IllegalArgumentException("Workflow terminated! Expected 5 files, but found ${fileList.size()} at ${params.bwa_index}")
		}
		//println "Valid BWA index files detected. Proceeding: ${fileList}"
		}
	}
	
def verify_minimap2_database () {
	println "Verifying minimap database"
	if (!params.minimap2_index) {
		throw new IllegalArgumentException("Please specify path to minimap2 database. Parameter minimap2_index is not specified or is null.")
		}
	ch_mm2_db = Channel.fromPath("${params.minimap2_index}.{mmi}").toList()
	ch_mm2_db.then { fileList ->
		if (fileList == null || fileList.isEmpty()) {
			throw new IllegalArgumentException("Workflow terminated! No valid files found for the specified minimap2 database path: ${params.minimap2_index}")
		}
		if (fileList.size() != 1) {
			throw new IllegalArgumentException("Workflow terminated! Expected 1 file, but found ${fileList.size()} at ${params.minimap2_index}")
		}
		}
	}
*/

// Process to verify BWA database
process verifyBwaDatabase {
    input:
    val bwaIndexPath

    output:
    val bwaDbValidated

    script:
    """
    // Check if the required BWA database files are available
    def filesExist = file("${bwaIndexPath}.{amb,ann,bwt,pac,sa}").exists()
    
    if (!filesExist) {
        throw new IllegalArgumentException("No valid BWA files found at: ${bwaIndexPath}")
    }
    
    // If valid, return success
    bwaDbValidated = true
    """
}

// Process to verify Minimap2 database
process verifyMinimap2Database {
    input:
    val minimap2IndexPath

    output:
    val mm2DbValidated

    script:
    """
    // Check if the required Minimap2 database file is available
    def fileExists = file("${minimap2IndexPath}.mmi").exists()
    
    if (!fileExists) {
        throw new IllegalArgumentException("No valid Minimap2 file found at: ${minimap2IndexPath}")
    }
    
    // If valid, return success
    mm2DbValidated = true
    """
}

/*
workflow CHECK_DATABASES {
    take:
    ch_sampledict
    
    main:
    // Convert the sample dictionary to a list and extract unique mappers
    ch_sampledict.toList().then { list ->
        println "Collected data: ${list}"

        // Extract unique mappers
        def mappers = list.collect { it.mapper }.unique()
        println "Unique mappers to be used: ${mappers}"

        // Conditional verification for databases
        if (mappers.contains('minimap2')) {
            println "Verifying minimap database"
            verifyMinimap2Database(params.minimap2_index)
        }
        
        if (mappers.contains('bwa')) {
            println "Verifying BWA database"
            verifyBwaDatabase(params.bwa_index)
        }
    }
}
*/

workflow CHECK_DATABASES {
    take:
    ch_sampledict  
    
    main:
		if (!params.minimap2_index) {
			throw new IllegalArgumentException("Please specify path to minimap2 database. Parameter minimap2_index is not specified or is null.")
			}
		
		
		ch_sampledict.toList().then { list ->
		println "Collected data: ${list}"
		
		
		// Extract unique mappers
		def mappers = list.collect { it.mapper }.unique()
		println "Unique mappers to be used: ${mappers}"
		
		
		// Verify databases
		if (mappers.contains('minimap2')) {
			println "Verifying minimap database"
			
			//verify_minimap2_database ()
		}
		
		if (mappers.contains('bwa')) {
			println "Verifying bwa database"
			//verify_bwa_database ()
		}
		
		}
		 	
}


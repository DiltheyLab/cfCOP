nextflow.enable.dsl = 2

include { BWA_MEM_NANO } from '../modules/bwa_mem'
include { BWA_MEM_ILLU } from '../modules/bwa_mem'

include { MINIMAP2_NANO } from '../modules/minimap2'
include { MINIMAP2_ILLU } from '../modules/minimap2'

include { SAMTOOLS_SORT_BYNAME } from '../modules/samtools_sort'
    
process SKIP_SORTING {
    input: 
    val sample
    
    output:
    tuple val(sample), path("${sample.id}_*byname.sam")
    
    when:
        sample.sorted_sam
    
    script:
    if (sample.mapper == "minimap2") {
    """
    ln -s "${sample.sorted_sam}" "${sample.id}_mm2_byname.sam"
    """
    }
    else if (sample.mapper == "bwa") {
    """
    ln -s "${sample.sorted_sam}" "${sample.id}_abwa_byname.sam"
    """
    }
    
}

    
workflow MAP {
    
    take:
        ch_sampledict
    
    main:
		
		ch_bwa_nano = BWA_MEM_NANO(ch_sampledict)
        ch_bwa_illu = BWA_MEM_ILLU(ch_sampledict)
        ch_mm2_nano = MINIMAP2_NANO(ch_sampledict)
        ch_mm2_illu = MINIMAP2_ILLU(ch_sampledict)
       
        ch_sam = ch_bwa_nano.concat(ch_bwa_illu,ch_mm2_nano,ch_mm2_illu)
        ch_sam_sorted = SAMTOOLS_SORT_BYNAME(ch_sam)
        ch_skipped_sorting = SKIP_SORTING(ch_sampledict)
        ch_sam_byname = ch_sam_sorted.concat(ch_skipped_sorting)
		
    emit:
        ch_sam_byname
}

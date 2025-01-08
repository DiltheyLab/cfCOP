nextflow.enable.dsl = 2

include { MINIMAP2_NANO } from '../modules/minimap2'
include { MINIMAP2_ILLU } from '../modules/minimap2'
include { SAMTOOLS_SORT_BYNAME } from '../modules/samtools_sort'

workflow MAP_MINIMAP2 {
    
    take:
        ch_sampledict
    
    main:
       ch_sam_nano = MINIMAP2_NANO(ch_sampledict)
       ch_sam_illu = MINIMAP2_ILLU(ch_sampledict)
       ch_sam = ch_sam_illu.concat(ch_sam_nano)
       ch_sam_sorted = SAMTOOLS_SORT_BYNAME(ch_sam)
       
   emit:
       ch_sam_sorted
}

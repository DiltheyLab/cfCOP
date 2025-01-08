nextflow.enable.dsl = 2

include { BWA_MEM_NANO } from '../modules/bwa_mem'
include { BWA_MEM_ILLU } from '../modules/bwa_mem'
include { SAMTOOLS_SORT_BYNAME } from '../modules/samtools_sort'

workflow MAP_BWA {
    
    take:
        ch_sampledict
    
    main:
       ch_sam_nano = BWA_MEM_NANO(ch_sampledict)
       ch_sam_illu = BWA_MEM_ILLU(ch_sampledict)
       ch_sam = ch_sam_illu.concat(ch_sam_nano)
       ch_sam_sorted = SAMTOOLS_SORT_BYNAME(ch_sam)
    
    emit:
       ch_sam_sorted
}

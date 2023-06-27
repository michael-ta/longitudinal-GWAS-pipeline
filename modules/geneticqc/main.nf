/* 
 * Tool log
 */
log.info """\
 LONG-GWAS - GENETICQC  P I P E L I N E
 ======================================
 chunk size for genetic data split and qc : ${params.chunk_size}
 assembly        : ${params.assembly}
 outdir       : ${PWD}/files/longGWAS_pipeline/results/cache/${params.dataset}
 """


/* 
 * Import modules
 */
include { GENETICQC } from './qc.nf'
include { MERGER_SPLITS; MERGER_CHRS } from './merge.nf'
//include { MERGER_CHRS } from './merge.nf'

workflow DOQC {
    take:
        input_check_ch
        cache
    
    main:
        input_check_ch
            .join(cache, remainder: true)
            .filter{ fSimple, fOrig, fCache -> fCache == null }
            .map{ fSimple, fOrig, fCache -> tuple( fSimple, fOrig ) }
            .set { chrvcf }

        chrvcf
            .count()
            .map{ val -> (1..val) }
            .flatten()
            .set{ input_idx_list }

        input_idx_list
            .merge(chrvcf)
            .branch {
                split_tasks_1: it[0] <= 5 
                                return  it
                split_tasks_2: it[0] < 10 && it[0] > 5 
                                return  it
                split_tasks_3: it[0] < 14 && it[0] >= 10 
                                return  it
                split_tasks_4: it[0] < 18 && it[0] >= 14 
                                return  it
                split_tasks_5: true
            }.set { input_p1 }

        input_p1.split_tasks_1
            .map{ fIdx, fSimple, fOrig -> fOrig }
            .splitText(by: params.chunk_size, file: true, compress: true)
            .map{ fn -> tuple( fn.getSimpleName(), fn) }
            .set{ input_split_ch_1 }

        input_p1.split_tasks_2
              .map{ fIdx, fSimple, fOrig -> fOrig }
              .splitText(by: params.chunk_size, file: true, compress: true)
              .map{ fn -> tuple( fn.getSimpleName(), fn) }
              .set{ input_split_ch_2 }

        input_p1.split_tasks_3
              .map{ fIdx, fSimple, fOrig -> fOrig }
              .splitText(by: params.chunk_size, file: true, compress: true)
              .map{ fn -> tuple( fn.getSimpleName(), fn) }
              .set{ input_split_ch_3 }

        input_p1.split_tasks_4
              .map{ fIdx, fSimple, fOrig -> fOrig }
              .splitText(by: params.chunk_size, file: true, compress: true)
              .map{ fn -> tuple( fn.getSimpleName(), fn) }
              .set{ input_split_ch_4 }

        input_p1.split_tasks_5
              .map{ fIdx, fSimple, fOrig -> fOrig }
              .splitText(by: params.chunk_size, file: true, compress: true)
              .map{ fn -> tuple( fn.getSimpleName(), fn) }
              .set{ input_split_ch_5 }

        input_split_ch_1
            .mix(input_split_ch_2, 
                 input_split_ch_3,
                 input_split_ch_4,
                 input_split_ch_5)
            .set{ input_chunks_ch }
        
        // input_split_ch_1
        //     .set{ input_chunks_ch }

        chrvcf
            .cross(input_chunks_ch)
            .flatten()
            .collate(4, false)
            .map{ fSimple1, fOrig, fSimple2, fSplit -> 
                tuple( fSimple1, file(fOrig), fSplit) }
            .set{ input_p1_run_ch }


        // RUN THE GENOTYPE DATA QC
        GENETICQC( input_p1_run_ch )
        GENETICQC.out.snpchunks_names.collectFile( newLine: true ) { vSimple, prefix ->
                                                                    ["${vSimple}.mergelist.txt", prefix] }
                                    .set{ chunknames }


        // DO ALL THE MERGING AFTER SPLITTING TO HANDLE BIG DATA MEMORY HUNTING
        MERGER_SPLITS(chunknames, GENETICQC.out.snpchunks_qc.collect())
        MERGER_SPLITS.out
            .collect()
            .flatten()
            .map{ fn -> tuple(fn.getSimpleName(), fn) }
            .concat( cache )
            .set{ chrsqced }

        chrsqced
            .collectFile() { vSimple, f ->
                [ "allchr.mergelist.txt", f.getBaseName() + '\n' ] }
            .set{ list_files_merge }
        chrsqced
            .map{ vSimple, f -> file(f) }
            .set{ chrfiles }

        MERGER_CHRS(list_files_merge  , chrfiles.collect())


        // Run the GWAS QC step
        

    emit:
        //MERGER.out | collect | flatten | map{ fn -> tuple(fn.getSimpleName(), fn) }
        //chrsqced | collect
        MERGER_CHRS.out | collect
}
/* 
 * Import modules
 */
include { GENETICQC }                           from '../modules/geneticqc/qc.nf'
include { MERGER_SPLITS; MERGER_CHRS }          from '../modules/geneticqc/merge.nf'
include { GWASQC }                              from '../modules/gwasqc/main.nf'


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
        
        chrvcf
            .cross(input_chunks_ch)
            .flatten()
            .collate(4, false)
            .map{ fSimple1, fOrig, fSimple2, fSplit -> 
                tuple( fSimple1, file(fOrig), fSplit) }
            .set{ input_p1_run_ch }


        // RUN THE GENOTYPE DATA QC
        GENETICQC( input_p1_run_ch )

        GENETICQC.out.snpchunks_names
            .collectFile( newLine: true ) 
                            { vSimple, prefix -> ["${vSimple}.mergelist.txt", prefix] }
            .set{ chunknames }

        // DO ALL THE MERGING AFTER SPLITTING TO HANDLE BIG DATA MEMORY HUNTING
        //MERGER_SPLITS(chunknames, GENETICQC.out.snpchunks_qc.collect())
        MERGER_SPLITS(chunknames, GENETICQC.out.snpchunks_merge.collect())
        
        MERGER_SPLITS.out
            .collect()
            .flatten()
            .map{ fn -> tuple(fn.getSimpleName(), fn) }
            .concat( cache )
            .set{ chrsqced }


        // Map some channels containing QC data that we will emit to use in other workflow        
        chrsqced
            .groupTuple(by: 0)
            .flatten()
            .collate(5)
            .set { gallop_plink_input }
        //gallop_plink_input.count().view()

        chrsqced
            .groupTuple(by: 0)
            .flatten()
            .filter( ~/.*pvar/ )
            .map { it -> tuple( it.getSimpleName(), it ) }
            .set{ gallopcph_chunks }
        //gallopcph_chunks.count().view()

        // Use chrsqced channel data to merge genotype data splits
        chrsqced
            .collectFile() { vSimple, f ->
                [ "allchr.mergelist.txt", f.getBaseName() + '\n' ] }
            .set{ list_files_merge }
        chrsqced
            .map{ vSimple, f -> file(f) }
            .set{ chrfiles }

        MERGER_CHRS(list_files_merge  , chrfiles.collect())
        MERGER_CHRS.out
            .flatten()
            .filter{ fName -> ["pgen", "pvar", "psam"].contains( fName.getExtension() ) }
            .collect()
            .set{  input_compute_pca }

        // Run the GWAS QC
        GWASQC(MERGER_CHRS.out)

        
    emit:
        chrchunks           = input_compute_pca  // Keep chr chunks we will use to compute PCA
        samplesqc           = GWASQC.out.gwasqc_h5_file 
        gallopcph_input     = gallopcph_chunks // this is the input pvar files we need to get ready for GALLOP and CPH
        gallopplink         = gallop_plink_input
}
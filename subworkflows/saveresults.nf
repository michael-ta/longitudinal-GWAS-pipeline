/*
 * Import modules
 */
include { SAVEGWAS }        from '../modules/resultsout/gwasfiles.nf'
include { MANHATTAN }       from '../modules/resultsout/manhattan.nf'


workflow SAVE_RESULTS {
    take:
        gwas

    main:

        SAVEGWAS(gwas)
        
        if ( params.mh_plot ) {
            MANHATTAN(SAVEGWAS.out.collect())
        }
}
/*
 * Import modules
 */
include { SAVEGWAS }        from '../modules/resultsout/gwasfiles.nf'
include { MANHATTAN }       from '../modules/resultsout/manhattan.nf'


workflow SAVE_RESULTS {
    take:
        gwas
        model
    main:
        SAVEGWAS(gwas, model)
        if ( params.mh_plot ) {
            MANHATTAN(SAVEGWAS.out.res_all.collect(), model)
        }
}
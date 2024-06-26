/*
 * Import modules
 */
include { GWASCPH }         from '../modules/gwasrun/cph.nf'
include { GWASGALLOP }      from '../modules/gwasrun/gallop.nf'
include { GWASGLM  }        from '../modules/gwasrun/glm.nf'


workflow GWAS_RUN {
    take:
        phenfile
        phennames
        gwaschunks
        glm_slist

    main:
        //phenoname_list = phennames.split(',')
        if ( params.longitudinal_flag) {
            GWASGALLOP(gwaschunks, phenfile, phennames)
            GWASRES = GWASGALLOP.out
        }
        else if ( params.survival_flag ) {
            GWASCPH(gwaschunks, phenfile, phennames)
            GWASRES = GWASCPH.out
        } else {
            GWASGLM(gwaschunks, glm_slist, phennames)
            //GWASGLM(gwaschunks, glm_slist)
            GWASRES = GWASGLM.out
        }

        GWASRES
            .groupTuple(sort: true)
            //.collect()
            .set { GROUP_RESULTS }
       //GWASRES
       //     .groupTuple()
       //     .map{ it[1] }
       //     .flatten()
       //    .collectFile(keepHeader: true)
       //    .set { RESULTS }

    emit:
        //RESULTS = RESULTS
        RESULTS = GROUP_RESULTS
}
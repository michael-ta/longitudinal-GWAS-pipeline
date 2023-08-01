/*
 * Import modules
 */
include { GWASCPH }         from '../modules/gwasrun/cph.nf'
include { GWASGALLOP }      from '../modules/gwasrun/gallop.nf'
include { GWASGLM  }        from '../modules/gwasrun/glm.nf'


workflow GWAS_RUN {
    take:
        phenos
        gwaschunks
        glm_slist

    main:
        if ( params.longitudinal_flag) {
            GWASGALLOP(gwaschunks, phenos)
            GWASRES = GWASGALLOP.out
        }
        else if ( params.survival_flag ) {
            GWASCPH(gwaschunks, phenos)
            GWASRES = GWASCPH.out
        } else {
            GWASGLM(gwaschunks, glm_slist)
            GWASRES = GWASGLM.out
        }

        GWASRES
            .groupTuple(sort: true)
            .collect()
            .set { GROUP_RESULTS }
       //GROUP_RESULTS.view()
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
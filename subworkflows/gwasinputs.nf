/*
 * Import modules
 */
include { GETPHENOS }                   from '../modules/gwasprep/outliers_exclude.nf'
include { REMOVEOUTLIERS }              from '../modules/gwasprep/outliers_exclude.nf'
include { COMPUTE_PCA; MERGE_PCA  }     from '../modules/gwasprep/covars.nf'
include { GALLOPCOX_INPUT }             from '../modules/gwasprep/gallopcph_in.nf'
include { RAWFILE_EXPORT }              from '../modules/gwasprep/raw.nf'
include { EXPORT_PLINK }                from '../modules/gwasprep/raw.nf'


workflow GWASDATA_PREP {
    take:
        covarfile
        qcworkflow

    main:
        PHENOS_FILE = GETPHENOS(covarfile)
        ALLPHENOS = GETPHENOS.out.allphenos.readLines()
        REMOVEOUTLIERS(qcworkflow.out.samplesqc, params.covarfile, ALLPHENOS)
        COMPUTE_PCA(REMOVEOUTLIERS.out.flatten(), qcworkflow.out.chrchunks)
        MERGE_PCA(COMPUTE_PCA.out.eigenvec)

        // If we are running survival or longitudinal, we format differently compared to GLM
        if  (params.longitudinal_flag | params.survival_flag) {
            GALLOPCOX_INPUT(qcworkflow.out.gallopcph_input)

            qcworkflow.out.gallopplink
                .cross(GALLOPCOX_INPUT.out.splitText(file: true))
                .flatten()
                .collate(7)
                .map{ it -> tuple(it[0], it[1], it[2], it[3], it[4], it[6]) }
                .set{ GALLOPCPHCHUNKS }

            RAWFILE_EXPORT(GALLOPCPHCHUNKS, MERGE_PCA.out)
            CHUNKS = RAWFILE_EXPORT.out.gwas_rawfile
            PLINK_SAMPLE_LIST = Channel.empty()

        } else {
            EXPORT_PLINK(MERGE_PCA.out.flatten(), params.phenofile)
            PLINK_SAMPLE_LIST = EXPORT_PLINK.out
            qcworkflow.out.gallopplink
                .set{ CHUNKS }
        }

    emit:
        CHUNKS          = CHUNKS
        PLINK_SLIST     = PLINK_SAMPLE_LIST
}
#!/usr/bin/env python

import pandas as pd
import time
import subprocess
import os
import argparse


def get_args():
    parser = argparse.ArgumentParser(description = "Gets phenotype and covar files, and remove some NAs for plink not to break")
    parser.add_argument("--pheno_covar", dest="pheno_covar", help="This is the unprocessed pheno and covar data", default=False)    
    parser.add_argument("--phenname", dest="phenname", help="This is the name of the single phenotype to run", default=False)
    parser.add_argument("--covname", dest="covname", help="Covariates to get")
    return parser.parse_args()


def phenocovar_glm(df, phenoname, covarname):
    # Read and generate covar and pheno files
    all_data = pd.read_csv(df, sep="\t", engine='c')
    all_data_filt = all_data.dropna(subset=[phenoname])
    pheno_file = all_data_filt.loc[:, ["#FID", "IID"] + [phenoname]]
    covar_file = all_data_filt.loc[:, ["#FID", "IID"] + covarname.split(' ') ]
    
    #Save the data
    pheno_file.to_csv("pheno.tsv", sep="\t", index=False)
    covar_file.to_csv("covar.tsv", sep="\t", index=False)

    print("covar and pheno files successfully generated")


if __name__ == "__main__":
    args = get_args()
    print(args)
    phenocovar_glm(df = args.pheno_covar, 
                    phenoname = args.phenname, 
                    covarname = args.covname)

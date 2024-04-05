#!/usr/bin/env python3

import pandas as pd
import os
import subprocess
import shutil
import numpy as np
import argparse
import plotly.express as px
import plotly
import plotly.io
plotly.io.kaleido.scope.chromium_args = tuple([arg for arg in plotly.io.kaleido.scope.chromium_args if arg != "--disable-dev-shm-usage"])

from qc import shell_do, het_prune, callrate_prune, merge_genos, get_outlier_ranges, ancestry_prune, king_prune, report_outliers, plot_3d

parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--ref', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format reference genotype file, everything before the *.bed/bim/fam.')
parser.add_argument('--ref_labels', type=str, default='nope', help='tab-separated plink-style IDs with ancestry label (FID  IID label) with no header')
parser.add_argument('--pop', type=str, default='nope', help='Population for analysis. right now, supports EUR and SAS')
parser.add_argument('--out', type=str, default='nope', help='Prefix for output (including path)')

args = parser.parse_args()

geno_path = args.geno
ref_path = args.ref
ref_labels = args.ref_labels
out_path = args.out
pop = args.pop

# make names for each step
callrate_out = f'{geno_path}_callrate'
het_out = f'{callrate_out}_het'
ancestry_out = f'{het_out}_ancestry'
king_out = f'{ancestry_out}_king'


# run steps
callrate = callrate_prune(geno_path, callrate_out, mind=0.05)

het = het_prune(callrate_out, het_out)
make_bed_cmd = f'plink2 --pfile {het_out} --make-bed --out {het_out}'
shell_do(make_bed_cmd)

ancestry = ancestry_prune(het_out, ref_path, ref_labels, ancestry_out, target_label=pop)

# make plots for pcs
total_pcs = pd.read_csv(f"{ancestry['output']['plink_out']}.pca.eigenvec", sep='\s+', header=None, names=['FID','IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'], dtype={'FID':str,'IID':str})
labels = pd.read_csv(ref_labels,sep='\t', header=None, names=['FID','IID','label'])
pcs_merge = total_pcs.merge(labels, on=['FID','IID'], how='left')
pcs_merge.label.fillna('new', inplace=True)

pc1_range = [pcs_merge.PC1.min(), pcs_merge.PC1.max()]
pc2_range = [pcs_merge.PC2.min(), pcs_merge.PC2.max()]
pc3_range = [pcs_merge.PC3.min(), pcs_merge.PC3.max()]

# plot reference panel + new sample pcs
plot_3d(pcs_merge, color='label', symbol=None, plot_out=f'{out_path}_total_pcs_plot', x='PC1', y='PC2', z='PC3', title='Reference Panel + New Sample PCs', x_range=pc1_range, y_range=pc2_range, z_range=pc3_range)

new_pcs = pcs_merge.loc[pcs_merge.label == 'new']
keep_samples = pd.read_csv(ancestry['output']['keep_samples'], sep='\t', dtype={'FID':str,'IID':str})
keep_samples.loc[:,'keep'] = True
new_pcs_keep = new_pcs.merge(keep_samples, how='left', on=['FID','IID'])
new_pcs_keep.loc[:,'keep'] = np.where(new_pcs_keep.keep == True, new_pcs_keep.keep, False)

# plot pcs for only new samples
plot_3d(new_pcs_keep, color='keep', symbol=None, plot_out=f'{out_path}_new_sample_pcs_plot', x='PC1', y='PC2', z='PC3', title='New Sample PCs', x_range=pc1_range, y_range=pc2_range, z_range=pc3_range)

eigenvals = pd.read_csv(f"{ancestry['output']['plink_out']}.pca.eigenval", header=None, names=['eigenval'])
eigenvals.loc[:,'PC'] = [f'PC{x+1}' for x in range(len(eigenvals.eigenval))]

fig = px.line(eigenvals, x='PC', y='eigenval')
fig.write_image(f'{out_path}_scree_plot.png', width=1980, height=1080)
fig.write_html(f'{out_path}_scree_plot.html')


make_pgen_cmd = f'plink2 --bfile {ancestry_out} --make-pgen --out {ancestry_out}'

shell_do(make_pgen_cmd)

king = king_prune(ancestry_out, king_out, king_cutoff=0.0884)


# put together list of outliers
steps = [callrate, het, ancestry, king]

outliers_path = f'{geno_path}.outliers'
outliers = report_outliers(steps, outliers_path)


# now add relevant dfs to hdf
ancestry_keep = pd.read_csv(ancestry['output']['keep_samples'], sep='\t')


# generate final pcs
pcs_out = f'{king_out}_pca'
pca_cmd = f'plink2 --pfile {king_out} --pca --out {pcs_out}'
shell_do(pca_cmd)
pcs = pd.read_csv(f'{pcs_out}.eigenvec', sep='\s+')

king_table_out = f'{king_out}_table'
king_table_cmd = f'plink2 --pfile {king_out} --make-king-table --out {king_table_out}'
shell_do(king_table_cmd)
kin = pd.read_csv(f'{king_table_out}.kin0', sep='\s+')

# add to hdf
outliers['outliers_df'].to_hdf(f'{out_path}.h5', key='outliers', mode='w')
ancestry_keep.to_hdf(f'{out_path}.h5', key='ancestry_keep')
pcs.to_hdf(f'{out_path}.h5', key='pcs')
kin.to_hdf(f'{out_path}.h5', key='kin')
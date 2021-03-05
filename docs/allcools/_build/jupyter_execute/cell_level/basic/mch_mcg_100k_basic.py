# Basic Cell Clustering Using 100Kb Bins

## Content

Here we go through the basic steps to perform cell clustering using genome non-overlapping 100Kb bins as features. We start from raw methylation data stored in MCDS format. It can be used to quickly evaluate get an idea on cell-type composition in a single-cell methylome dataset (e.g., the dataset from a single experiment).

## Input
- MCDS files
- Cell metadata

## Output
- Cell-by-100kb-bin AnnData with embedding coordinates and cluster labels.

## Import

%load_ext autoreload

%autoreload 2

import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

from ALLCools import MCDS
from ALLCools.clustering import tsne, significant_pc_test
from ALLCools.plot import *

## Parameters

# change this to the path to your metadata
metadata_path = '../../data/MOp/MOp.CellMetadata.csv.gz'

# Basic filtering parameters
mapping_rate_cutoff = 0.5
mapping_rate_col_name = 'MappingRate'  # Name may change
final_reads_cutoff = 500000
final_reads_col_name = 'FinalmCReads'  # Name may change
mccc_cutoff = 0.03
mccc_col_name = 'mCCCFrac'  # Name may change
mch_cutoff = 0.2
mch_col_name = 'mCHFrac'  # Name may change
mcg_cutoff = 0.5
mcg_col_name = 'mCGFrac'  # Name may change

# change this to the paths to your MCDS files, 
# ALLCools.MCDS can handle multiple MCDS files automatically
mcds_path_list = [
    '../../data/MOp/3C-171206.mcds'
]

# Dimension name used to do clustering
# This corresponding to AnnData .obs and .var
obs_dim = 'cell'  # observation
var_dim = 'chrom100k'  # feature

# feature cov cutoffs
min_cov = 500
max_cov = 3000

# Regions to remove during the clustering analysis
# change this to the path to ENCODE blacklist.
# The ENCODE blacklist can be download from https://github.com/Boyle-Lab/Blacklist/
black_list_path = '../../data/genome/mm10-blacklist.v2.bed.gz'
black_list_fraction = 0.2
exclude_chromosome = ['chrM', 'chrY']

# load to memory or not
load = True

# HVF
mch_pattern = 'CHN'
mcg_pattern = 'CGN'
n_top_feature = 20000

# PC cutoff
pc_cutoff = 0.1

# KNN
knn = -1  # -1 means auto determine

# Leiden
resolution = 1

## Load Cell Metadata

metadata = pd.read_csv(metadata_path, index_col=0)
print(f'Metadata of {metadata.shape[0]} cells')
metadata.head()

## Filter Cells

judge = (metadata[mapping_rate_col_name] > mapping_rate_cutoff) & \
        (metadata[final_reads_col_name] > final_reads_cutoff) & \
        (metadata[mccc_col_name] < mccc_cutoff) & \
        (metadata[mch_col_name] < mch_cutoff) & \
        (metadata[mcg_col_name] > mcg_cutoff)

metadata = metadata[judge].copy()
print(f'{metadata.shape[0]} cells passed filtering')

metadata.to_csv('CellMetadata.PassQC.csv.gz')

## Load MCDS

mcds = MCDS.open(
    mcds_path_list, 
    obs_dim='cell', 
    use_obs=metadata.index  # MCDS contains all cells, this will select cells that passed filtering 
)
total_feature = mcds.get_index(var_dim).size
mcds

# you can add the cell metadata into MCDS
mcds.add_cell_metadata(metadata)

## Filter Features

mcds.add_feature_cov_mean(var_dim=var_dim)

We saw three parts here, from coverage low to high, they are
1. Low coverage regions
2. chrX regions, because this dataset from male mouse brain
3. Other autosomal regions

# filter by coverage - based on the distribution above
mcds = mcds.filter_feature_by_cov_mean(
    var_dim=var_dim,
    min_cov=min_cov,  # minimum coverage
    max_cov=max_cov  # Maximum coverage
)

# remove blacklist regions
mcds = mcds.remove_black_list_region(
    var_dim,
    black_list_path,
    f=black_list_fraction  # Features having overlap > f with any black list region will be removed.
)

# remove chromosomes
mcds = mcds.remove_chromosome(var_dim, exclude_chromosome)

## Calculate Feature mC Fractions

mcds.add_mc_frac(
    var_dim=var_dim, 
    normalize_per_cell=True,  # after calculating mC frac, per cell normalize the matrix
    clip_norm_value=10  # clip outlier values above 10 to 10
)

# load only the mC fraction matrix into memory so following steps is faster
# Only load into memory when you memory size is enough to handle your dataset
if load and (mcds.get_index(obs_dim).size < 20000):
    mcds[f'{var_dim}_da_frac'].load()

The RuntimeWarning is expected (due to cov == 0). You can ignore it.

## Select Highly Variable Features (HVF)

### mCH HVF

mch_hvf = mcds.calculate_hvf_svr(var_dim=var_dim,
                                 mc_type=mch_pattern,
                                 n_top_feature=n_top_feature,
                                 plot=True)

### mCG HVF

mcg_hvf = mcds.calculate_hvf_svr(var_dim=var_dim,
                                 mc_type=mcg_pattern,
                                 n_top_feature=n_top_feature,
                                 plot=True)

## Get cell-by-feature mC fraction AnnData

mch_adata = mcds.get_adata(mc_type=mch_pattern,
                           var_dim=var_dim,
                           select_hvf=False)
mch_adata

mcg_adata = mcds.get_adata(mc_type=mcg_pattern, var_dim=var_dim, select_hvf=False)
mcg_adata

## PCA

### mCH PCA

sc.tl.pca(mch_adata)
ch_n_components = significant_pc_test(mch_adata)
fig, axes = plot_decomp_scatters(mch_adata,
                                 n_components=ch_n_components,
                                 hue=mch_col_name,
                                 hue_quantile=(0.25, 0.75),
                                 nrows=3,
                                 ncols=5)

### mCG PCA

sc.tl.pca(mcg_adata)
cg_n_components = significant_pc_test(mcg_adata)
fig, axes = plot_decomp_scatters(mcg_adata,
                                 n_components=cg_n_components,
                                 hue=mch_col_name,
                                 hue_quantile=(0.25, 0.75),
                                 nrows=3,
                                 ncols=5)

### Concatenate PCs

ch_pcs = mch_adata.obsm['X_pca'][:, :ch_n_components]
cg_pcs = mcg_adata.obsm['X_pca'][:, :cg_n_components]

# scale the PCs so CH and CG PCs has the same total var
cg_pcs = cg_pcs / cg_pcs.std()
ch_pcs = ch_pcs / ch_pcs.std()

# total_pcs
total_pcs = np.hstack([ch_pcs, cg_pcs])

# make a copy of adata, add new pcs
# this is suboptimal, will change this when adata can combine layer and X in the future
adata = mch_adata.copy()
adata.obsm['X_pca'] = total_pcs
del adata.uns['pca']
del adata.varm['PCs']

## Clustering

### Calculate Nearest Neightbors

if knn == -1:
    knn = max(15, int(np.log2(adata.shape[0])*2))
sc.pp.neighbors(adata, n_neighbors=knn)

### Leiden Clustering

sc.tl.leiden(adata, resolution=resolution)

## Manifold learning

def dump_embedding(adata, name, n_dim=2):
    # put manifold coordinates into adata.obs
    for i in range(n_dim):
        adata.obs[f'{name}_{i}'] = adata.obsm[f'X_{name}'][:, i]
    return adata

### tSNE

tsne(adata,
     obsm='X_pca',
     metric='euclidean',
     exaggeration=-1,  # auto determined
     perplexity=30,
     n_jobs=-1)
adata = dump_embedding(adata, 'tsne')

fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = categorical_scatter(data=adata.obs,
                        ax=ax,
                        coord_base='tsne',
                        hue='leiden',
                        text_anno='leiden',
                        show_legend=True)

### UMAP

sc.tl.umap(adata)
adata = dump_embedding(adata, 'umap')

fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = categorical_scatter(data=adata.obs,
                        ax=ax,
                        coord_base='umap',
                        hue='leiden',
                        text_anno='leiden',
                        show_legend=True)

### Interactive plot

interactive_scatter(data=adata.obs, hue='leiden', coord_base='umap')

## Save Results

adata.write_h5ad('adata.chrom100k-clustering.h5ad')
adata

adata.obs.to_csv('clustering_results.csv.gz')
adata.obs.head()

## Sanity test

# This test dataset come from Liu et al. 2021 Nature, so we already annotated the cell types
# For new datasets, see following notebooks about identifying cluster markers and annotate clusters
if 'CellTypeAnno' in adata.obs:
    adata.obs['CellTypeAnno'] = adata.obs['CellTypeAnno'].fillna('Outlier')
    
    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
    _ = categorical_scatter(data=adata.obs,
                            ax=ax,
                            coord_base='umap',
                            hue='CellTypeAnno',
                            text_anno='CellTypeAnno',
                            palette='tab20',
                            show_legend=True)

You may notice that here is an outlier population near the IT-L4, IT-L5, IT-L6, which is likely correspond to potential doublets. We can identify doublets using the MethylScrublet notebook.


# Basic Cell Clustering Using mCG-5Kb Bins

## Content

Here we go through the basic steps to perform cell clustering using genome non-overlapping 5Kb bins as features. We start from hypo-methylation probability data stored in mC-AnnData format (.mcad to distinguish from normal .h5ad files). It can be used to quickly evaluate get an idea on cell-type composition in a single-cell methylome dataset (e.g., the dataset from a single experiment). Comparing to 100Kb bins clustering process, this clustering process is more suitable for samples with low mCH fraction (many non-brain tissues) and narrow methylation diversity (so smaller feature works better).

### Dataset used in this notebook
- Adult (age P56) male mouse brain puititary (PIT) snmC-seq2 data from Frederique et al. 2021 (REF)

## Input
- MCAD file
- Cell metadata

## Output
- Cell-by-5kb-bin AnnData (sparse matrix) with embedding coordinates and cluster labels.

## Import

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata
import scanpy as sc

from ALLCools.clustering import tsne, significant_pc_test, filter_regions, remove_black_list_region, tf_idf, binarize_matrix
from ALLCools.plot import *

## Parameters

metadata_path = '../../data/PIT/PIT.CellMetadata.csv.gz'
mcad_path = '../../data/PIT/PIT.mcad'

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
# cell metadata for this example is filtered already
print(f'{metadata.shape[0]} cells passed filtering')

## Load MCAD

mcad = anndata.read_h5ad(mcad_path)
# add cell metadata to mcad:
mcad.obs = pd.concat([mcad.obs, metadata.reindex(mcad.obs_names)], axis=1)
mcad

## Binarize

binarize_matrix(mcad, cutoff=0.95)

## Filter Features

filter_regions(mcad, hypo_cutoff=6)

remove_black_list_region(mcad, black_list_path='/home/hanliu/ref/blacklist/mm10-blacklist.v2.bed.gz')

mcad

## TF-IDF Transform and Dimension Reduction

# by default we save the results in adata.obsm['X_pca'] which is the scanpy defaults in many following functions
# But this matrix is not calculated by PCA
tf_idf(mcad, algorithm='arpack', obsm='X_pca')

# choose significant components
significant_pc_test(mcad, p_cutoff=pc_cutoff, update=True)

## Clustering

### Calculate Nearest Neighbors

if knn == -1:
    knn = max(15, int(np.log2(mcad.shape[0])*2))
sc.pp.neighbors(mcad, n_neighbors=knn)

## Leiden Clustering

sc.tl.leiden(mcad, resolution=resolution)

## Manifold learning

def dump_embedding(adata, name, n_dim=2):
    # put manifold coordinates into adata.obs
    for i in range(n_dim):
        adata.obs[f'{name}_{i}'] = adata.obsm[f'X_{name}'][:, i]
    return adata

### tSNE

tsne(mcad,
     obsm='X_pca',
     metric='euclidean',
     exaggeration=-1,  # auto determined
     perplexity=30,
     n_jobs=-1)
mcad = dump_embedding(mcad, 'tsne')

fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = categorical_scatter(data=mcad.obs,
                        ax=ax,
                        coord_base='tsne',
                        hue='leiden',
                        text_anno='leiden',
                        show_legend=True)

### UMAP

sc.tl.umap(mcad)
mcad = dump_embedding(mcad, 'umap')

fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = categorical_scatter(data=mcad.obs,
                        ax=ax,
                        coord_base='umap',
                        hue='leiden',
                        text_anno='leiden',
                        show_legend=True)

### Interactive plot

interactive_scatter(data=mcad.obs, hue='leiden', coord_base='umap')

## Save results

mcad.write_h5ad('PIT.mCG-5K-clustering.h5ad')
mcad

mcad.obs.to_csv('PIT.ClusteringResults.csv.gz')
mcad.obs.head()

## Sanity test

This test dataset come from Frederique et al. 2021 (REF), so we already annotated the cell types. For new datasets, see following notebooks about identifying cluster markers and annotate clusters

if 'CellTypeAnno' in mcad.obs:
    mcad.obs['CellTypeAnno'] = mcad.obs['CellTypeAnno'].fillna('Outlier')
    
    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
    _ = categorical_scatter(data=mcad.obs,
                            ax=ax,
                            coord_base='umap',
                            hue='CellTypeAnno',
                            text_anno='CellTypeAnno',
                            palette='tab20',
                            show_legend=True)


# Consensus Clustering

## Purpose
The purpose of this step is to run consensus clustering.

## Input
- HVF adata files.

## Output
- HVF adata file with cluster annotated.

## Import

import pathlib
import anndata
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from ALLCools.clustering import ConsensusClustering, Dendrogram
from ALLCools.plot import *

## Parameters

# clustering name
clustering_name = 'L1'

# input data
cell_meta_path = './CellMetadata.PassQC.csv.gz'
adata_path = './adata.with_coords.h5ad'
coord_base = 'tsne'

# ConsensusClustering
n_neighbors = 25
metric = 'euclidean'
min_cluster_size = 10
consensus_rate = 0.5
leiden_repeats = 500
leiden_resolution = 0.8
random_state = 0
n_jobs = 40
train_frac = 0.5
train_max_n = 500
max_iter = 20
outlier_cell_cutoff = 0.25
outlier_label = 'Outlier'

# Dendrogram via Multiscale Bootstrap Resampling
nboot = 10000
method_dist = 'correlation'
method_hclust = 'average'

plot_type = 'static'

## Load Data

cell_meta = pd.read_csv(cell_meta_path, index_col=0)
adata = anndata.read_h5ad(adata_path)

## Consensus Clustering

cc = ConsensusClustering(model=None,
                         n_neighbors=n_neighbors,
                         metric=metric,
                         min_cluster_size=min_cluster_size,
                         consensus_rate=consensus_rate,
                         repeats=leiden_repeats,
                         resolution=leiden_resolution,
                         random_state=random_state,
                         train_frac=train_frac,
                         train_max_n=train_max_n,
                         max_iter=max_iter,
                         n_jobs=n_jobs,
                         outlier_cell_cutoff=outlier_cell_cutoff,
                         outlier_label=outlier_label,
                         plot=True)

if 'X_pca' not in adata.obsm:
    raise KeyError(
        'X_pca do not exist in the adata file, run PCA first before clustering.'
    )
cc.fit_predict(adata.obsm['X_pca'])

## Plot

### Cluster Lables

adata.obs[clustering_name] = cc.label

fig, ax = plt.subplots(figsize=(4, 4), dpi=250)
_ = categorical_scatter(data=adata.obs,
                        ax=ax,
                        hue=clustering_name,
                        coord_base=coord_base,
                        palette='tab20',
                        text_anno=clustering_name,
                        show_legend=True)

### Final Prediction Probability

adata.obs[clustering_name + '_proba'] = cc.label_proba

fig, ax = plt.subplots(figsize=(4, 4), dpi=250)
_ = continuous_scatter(data=adata.obs,
                       ax=ax,
                       hue_norm=(0, 1),
                       hue=clustering_name + '_proba',
                       coord_base=coord_base)

### Prediction Probability Per Cluster

fig, ax = plt.subplots(figsize=(6, 3), dpi=300)

sns.violinplot(data=adata.obs,
               x=clustering_name,
               y=clustering_name + '_proba',
               scale='width',
               linewidth=0.5,
               cut=0,
               ax=ax)
ax.set(ylim=(0, 1), title='Prediction Probability Per Cluster')
ax.xaxis.set_tick_params(rotation=90)
ax.grid(linewidth=0.5, color='gray', linestyle='--')
sns.despine(ax=ax)

## Calculate Cluster Dendrogram

# using the cluster centroids in PC space to calculate dendrogram
pc_matrix = adata.obsm['X_pca']
pc_center = pd.DataFrame(pc_matrix, index=adata.obs_names).groupby(
    adata.obs[clustering_name]).median()
pc_center = pc_center[pc_center.index != outlier_label]
# Dendrogram take feature-by-sample dataframe as R pvclust function does.
pc_center = pc_center.T
idx_to_cluster = {i: c for i, c in enumerate(pc_center.index)}

pc_center.to_hdf('ClusterPCCenters.hdf', key='data')
pc_center.shape

dendro = Dendrogram(nboot=nboot,
                    method_dist=method_dist,
                    method_hclust=method_hclust,
                    n_jobs=n_jobs)
dendro.fit(pc_center, plot=True)

## Sanity Test

if 'CellTypeAnno' in cell_meta:
    fig, axes = plt.subplots(figsize=(8, 4), dpi=250, ncols=2)
    ax = axes[0]
    _ = categorical_scatter(data=adata.obs,
                            ax=ax,
                            hue=clustering_name,
                            coord_base=coord_base,
                            palette='tab20',
                            text_anno=clustering_name,
                            show_legend=False)
    ax = axes[1]
    adata.obs['CellTypeAnno'] = cell_meta['CellTypeAnno']
    _ = categorical_scatter(data=adata.obs.dropna(subset=['CellTypeAnno']),
                            ax=ax,
                            hue='CellTypeAnno',
                            coord_base=coord_base,
                            palette='tab20',
                            text_anno='CellTypeAnno',
                            show_legend=False)

## Save

cc.save('ConcensusClustering.model.lib')
dendro.save('Dendrogram.lib')
adata.write_h5ad(adata_path)

adata


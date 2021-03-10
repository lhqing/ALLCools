# Plot Cluster DMGs

import pandas as pd
import scanpy as sc
import anndata
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
import dask
from ALLCools.plot import *
from ALLCools.mcds.mcds import MCDS
import pathlib
import numpy as np
import warnings
import joblib
import networkx as nx

clustering_name = 'L1'
downsample = 30000
gene_fraction_dir = 'gene_frac/'
mc_type = 'CHN'
coord_base = 'tsne'
gene_meta_path = '/home/hanliu/ref/mouse/gencode/vm22/gencode.vM22.annotation.gene.flat.tsv.gz'
ncols = 5
nrows = 2
outlier_label = 'Outlier'
load_top_n = 50
adata_path = 'adata.coord_only.h5ad'

## Load

### Clustering results

adata = anndata.read_h5ad(adata_path)

if downsample and (adata.n_obs > downsample):
    use_cells = adata.obs.sample(downsample, random_state=0).index
    adata = adata[adata.obs_names.isin(use_cells), :].copy()
else:
    use_cells = adata.obs_names

## Cluster DMGs

cluster_dmgs = {}
with pd.HDFStore('ClusterRankedDMG.CHN.hdf') as hdf:
    for cluster in hdf.keys():
        cluster_dmgs[cluster[1:]] = hdf[cluster][:load_top_n]

# only load relevant genes
use_genes = set()
for k in cluster_dmgs.values():
    use_genes.update(set(k.index))

## Gene mC Fraction Data

gene_meta = pd.read_csv(gene_meta_path, sep='\t', index_col='gene_id')

gene_fraction_dir = pathlib.Path(gene_fraction_dir)
gene_meta = pd.read_csv(gene_fraction_dir / 'GeneMetadata.csv.gz', index_col=0)
gene_meta.index.name = 'gene_id'
gene_frac_da = xr.open_mfdataset(f'{gene_fraction_dir}/*_da_rate.nc',
                                 concat_dim='cell',
                                 combine='nested')[f'gene_da_rate']

# chose cell and genes to load
_cells = gene_frac_da.get_index('cell')
_cells = _cells[_cells.isin(use_cells)]
_genes = gene_frac_da.get_index('gene')
_genes = _genes[_genes.isin(use_genes)]
gene_frac_da = gene_frac_da.sel(cell=_cells, gene=_genes, mc_type=mc_type).load()
gene_meta = gene_meta.loc[_genes].copy()
gene_frac_da

## Plot

cell_meta = adata.obs.copy()

### Plot Function

def plot_genes(cluster, ncols, nrows, axes_size=3):
    ncols = max(2, ncols)
    nrows += 1
    
    # figure
    fig = plt.figure(figsize=(ncols * axes_size, nrows * axes_size), dpi=150)
    gs = fig.add_gridspec(nrows=nrows, ncols=ncols)
    
    # cluster axes
    ax = fig.add_subplot(gs[0, 0])
    categorical_scatter(data=cell_meta,
                        ax=ax,
                        coord_base=coord_base,
                        axis_format=None,
                        hue=clustering_name,
                        palette='tab20')
    ax.set_title('All Clusters')
    ax = fig.add_subplot(gs[0, 1])
    categorical_scatter(data=cell_meta,
                        ax=ax,
                        coord_base=coord_base,
                        hue=cell_meta[clustering_name] == cluster,
                        axis_format=None,
                        palette={
                            True: 'red',
                            False: 'lightgray'
                        })
    ax.set_title('This Cluster')
    
    # gene axes
    axes = []
    for row in range(1, nrows):
        for col in range(ncols):
            axes.append(fig.add_subplot(gs[row, col]))
    n_genes = len(axes)
    
    for ax, (gene, value) in zip(axes, cluster_dmgs[cluster][:n_genes].items()):
        if ax.is_first_col() and ax.is_last_row():
            axis = 'tiny'
        else:
            axis = None
        continuous_scatter(ax=ax,
                           data=cell_meta,
                           hue=gene_frac_da.sel(gene=gene).to_pandas(),
                           axis_format=axis,
                           hue_norm=(0.67, 1.5),
                           coord_base=coord_base)
        name = gene_meta.loc[gene, 'gene_name']
        ax.set_title(f'{name} ({value:.2f})')
    fig.suptitle(f'Top {cluster} Markers')
    return

### Per Cluster Plots

for cluster in sorted(cell_meta[clustering_name].unique()):
    if cluster != outlier_label:
        plot_genes(cluster, ncols, nrows, axes_size=3)


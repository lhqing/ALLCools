# Differential Methylated Genes - Pairwise

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
from ALLCools.clustering import PairwiseDMG, aggregate_pairwise_dmg
import pathlib
import numpy as np
from itertools import combinations
import warnings
from sklearn.metrics import roc_auc_score
from concurrent.futures import ProcessPoolExecutor, as_completed

gene_meta_path = '/home/hanliu/ref/mouse/gencode/vm22/gencode.vM22.annotation.gene.flat.tsv.gz'
chrom_to_remove = ['chrM']
adata_path = 'adata.coord_only.h5ad'
clustering_name = 'L1'

# change this to the path to your filtered metadata
metadata_path = 'CellMetadata.PassQC.csv.gz'

# change this to the paths to your MCDS files
gene_fraction_dir = 'gene_frac/'
obs_dim = 'cell'

# DMG
mc_type = 'CHN'
min_cluster_cell_number = 10
top_n = 30000
max_cell_per_group = 1000
chunk_size = 100
cpu = 10
random_state = 0
adj_p_cutoff = 1e-3
delta_rate_cutoff = 0.3
auroc_cutoff = 0.9
n_jobs = 30

## Load

### Gene Fraction Data

gene_fraction_dir = pathlib.Path(gene_fraction_dir)
gene_meta = pd.read_csv('GeneMetadata.csv.gz', index_col=0)
gene_meta.index.name = 'gene_id'

gene_frac_da = xr.open_mfdataset(f'{gene_fraction_dir}/*_da_rate.nc',
                                 concat_dim='cell',
                                 combine='nested')[f'gene_da_rate']
# standardize names
gene_frac_da

### Clustering Data

metadata = pd.read_csv(metadata_path, index_col=0)
total_cells = metadata.shape[0]
print(f'Metadata of {total_cells} cells')

adata = anndata.read_h5ad(adata_path)
adata

## Pairwise DMG

pwdmg = PairwiseDMG(max_cell_per_group=max_cell_per_group,
                    top_n=top_n,
                    adj_p_cutoff=adj_p_cutoff,
                    delta_rate_cutoff=delta_rate_cutoff,
                    auroc_cutoff=auroc_cutoff,
                    random_state=random_state,
                    n_jobs=n_jobs)

pwdmg.fit_predict(x=gene_frac_da.sel(mc_type=mc_type, cell=adata.obs_names), 
                  groups=adata.obs[clustering_name],
                  cleanup=True,
                  outlier='Outlier')

pwdmg.dmg_table.to_hdf(f'PairwiseDMG.{mc_type}.hdf', key='data')

## Cluster DMG

cluster_dmgs = aggregate_pairwise_dmg(pwdmg.dmg_table, adata, groupby=clustering_name)

with pd.HDFStore(f'ClusterRankedDMG.{mc_type}.hdf') as hdf:
    for cluster, dmgs in cluster_dmgs.items():
        hdf[cluster] = dmgs


# Calculate Highly Variable Features And Get mC Fraction AnnData

## Purpose
The purpose of this step is to select highly variable features (HVF) and generate cell-by-feature methylation fraction matrix for clustering. The highly variable features are selected by comparing feature's normalized dispersion among cells.

## Input
- Filtered cell metadata;
- MCDS files;
- Feature list from basic feature filtering

## Output
- cell-by-HVF methylation fraction matrix stored in AnnData format, e.g., mCH adata and mCG adata.

## Import

import pandas as pd
import xarray as xr
import dask
import ALLCools
from ALLCools import MCDS
from ALLCools.clustering import feature_enrichment

## Parameters

# If True, will load all data into memory.
# Computation will be much faster, but also very memory intensive, only use this for small number of cells (<10,000)
load = True

# change this to the path to your filtered metadata
metadata_path = 'CellMetadata.PassQC.csv.gz'

# change this to the paths to your MCDS files
mcds_path_list = [
    '../../../data/Brain/3C-171206.mcds',
    '../../../data/Brain/3C-171207.mcds',
    '../../../data/Brain/9H-190212.mcds',
    '../../../data/Brain/9H-190219.mcds',
]

# Feature list after basic filter
feature_path = 'FeatureList.BasicFilter.txt'

# Dimension name used to do clustering
obs_dim = 'cell'  # observation
var_dim = 'chrom100k'  # feature

# HVF method:
# SVR: regression based
# Bins: normalize dispersion per bin
hvf_method = 'SVR'
mch_pattern = 'CHN'
mcg_pattern = 'CGN'
n_top_feature = 20000

# Downsample cells
downsample = 20000

## Load Data

### Metadata

metadata = pd.read_csv(metadata_path, index_col=0)
total_cells = metadata.shape[0]
print(f'Metadata of {total_cells} cells')

metadata.head()

use_features = pd.read_csv(feature_path, header=None, index_col=0).index
use_features.name = var_dim

### MCDS

with dask.config.set(**{'array.slicing.split_large_chunks': False}):
    # still use all the cells to load MCDS
    total_mcds = MCDS.open(mcds_path_list,
                           obs_dim=obs_dim,
                           use_obs=metadata.index).sel({var_dim: use_features})

## Add mC Rate

total_mcds.add_mc_rate(var_dim=var_dim,
                       normalize_per_cell=True,
                       clip_norm_value=10)

total_mcds

### If downsample

if downsample and total_cells > downsample:
    # make a downsampled mcds
    print(f'Downsample cells to {downsample} to calculate HVF.')
    downsample_cell_ids = metadata.sample(downsample, random_state=0).index
    mcds = total_mcds.sel(
        {obs_dim: total_mcds.get_index(obs_dim).isin(downsample_cell_ids)})
else:
    mcds = total_mcds

if load and (mcds.get_index('cell').size <= 20000):
    # load the relevant data so the computation can be fater, watch out memory!
    mcds[f'{var_dim}_da_frac'].load()

The RuntimeWarning is expected (due to cov == 0). You can ignore it.

## Highly Variable Feature

### mCH

if hvf_method == 'SVR':
    # use SVR based method
    mch_hvf = mcds.calculate_hvf_svr(var_dim=var_dim,
                                     mc_type=mch_pattern,
                                     n_top_feature=n_top_feature,
                                     plot=True)
else:
    # use bin based method
    mch_hvf = mcds.calculate_hvf(var_dim=var_dim,
                                 mc_type=mch_pattern,
                                 min_mean=0,
                                 max_mean=5,
                                 n_top_feature=n_top_feature,
                                 bin_min_features=5,
                                 mean_binsize=0.05,
                                 cov_binsize=100)

### mCG

if hvf_method == 'SVR':
    # use SVR based method
    mcg_hvf = mcds.calculate_hvf_svr(var_dim=var_dim,
                                     mc_type=mcg_pattern,
                                     n_top_feature=n_top_feature,
                                     plot=True)
else:
    # use bin based method
    mcg_hvf = mcds.calculate_hvf(var_dim=var_dim,
                                 mc_type=mcg_pattern,
                                 min_mean=0,
                                 max_mean=5,
                                 n_top_feature=n_top_feature,
                                 bin_min_features=5,
                                 mean_binsize=0.02,
                                 cov_binsize=20)

## Save AnnData

mch_adata = mcds.get_adata(mc_type=mch_pattern,
                           var_dim=var_dim,
                           select_hvf=True)

mch_adata.write_h5ad(f'mCH.HVF.h5ad')

mch_adata

mcg_adata = mcds.get_adata(mc_type=mcg_pattern,
                           var_dim=var_dim,
                           select_hvf=True)

mcg_adata.write_h5ad(f'mCG.HVF.h5ad')

mcg_adata


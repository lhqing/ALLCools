# Feature Basic Filtering

## Purpose
Apply basic filters to remove these problematic features:
- Extremly low coverage or high coverage features
- ENCODE Blcaklist
- Some chromosomes (usually, chrY and chrM)

## Input
- Cell metadata (after basic cell filter)
- MCDS files

## Output
- FeatureList.BasicFilter.txt: List of feature ids passed all filters

## Import

import pandas as pd
import seaborn as sns
from ALLCools import MCDS

sns.set_context(context='notebook', font_scale=1.3)

## Parameters

# change this to the path to your filtered metadata
metadata_path = 'CellMetadata.PassQC.csv.gz'

# change this to the paths to your MCDS files
mcds_path_list = [
    '../../../data/Brain/3C-171206.mcds',
    '../../../data/Brain/3C-171207.mcds',
    '../../../data/Brain/9H-190212.mcds',
    '../../../data/Brain/9H-190219.mcds',
]

# Dimension name used to do clustering
obs_dim = 'cell'  # observation
var_dim = 'chrom100k'  # feature

min_cov = 500
max_cov = 3000

# change this to the path to ENCODE blacklist.
# The ENCODE blacklist can be download from https://github.com/Boyle-Lab/Blacklist/
black_list_path = '../../../data/genome/mm10-blacklist.v2.bed.gz'
f = 0.2

exclude_chromosome = ['chrM', 'chrY']

## Load Data

### Metadata

metadata = pd.read_csv(metadata_path, index_col=0)
total_cells = metadata.shape[0]
print(f'Metadata of {total_cells} cells')

metadata.head()

### MCDS

mcds = MCDS.open(mcds_path_list, obs_dim='cell', use_obs=metadata.index)
total_feature = mcds.get_index(var_dim).size

mcds

## Filter Features

### Filter by mean coverage

mcds.add_feature_cov_mean(var_dim=var_dim)

mcds = mcds.filter_feature_by_cov_mean(
    var_dim=var_dim,
    min_cov=min_cov,  # minimum coverage
    max_cov=max_cov  # Maximum coverage
)

### Filter by ENCODE Blacklist

mcds = mcds.remove_black_list_region(
    var_dim,
    black_list_path,
    f=f  # Features having overlap > f with any black list region will be removed.
)

### Remove chromosomes

mcds = mcds.remove_chromosome(var_dim, exclude_chromosome)

## Save Feature List

print(
    f'{mcds.get_index(var_dim).size} ({mcds.get_index(var_dim).size * 100 / total_feature:.1f}%) '
    f'{var_dim} remained after all the basic filter.')

with open('FeatureList.BasicFilter.txt', 'w') as f:
    for var in mcds.get_index(var_dim).astype(str):
        f.write(var + '\n')


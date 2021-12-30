import numpy as np
import yaml


def calculate_chunk_regions(n_features, dtype, chunk_size_gbs):
    # calculate regions per chunk
    chunk_mem_size = chunk_size_gbs * 1024 * 1024 * 1024
    dtype_size = np.dtype(dtype).itemsize
    n_regions = int(np.floor(chunk_mem_size / n_features / dtype_size / 10000) * 10000)
    n_regions = max(n_regions, 1000)
    n_regions = min(n_regions, 200000)
    return n_regions


def update_region_ds_config(output_dir, new_dataset_dim=None, change_region_dim=None):
    # update RegionDS default dimension
    with open(f'{output_dir}/.ALLCools', 'r') as f:
        config = yaml.load(f, yaml.SafeLoader)
    if change_region_dim is not None:
        config['region_dim'] = change_region_dim
    if new_dataset_dim is not None:
        config['ds_region_dim'].update(new_dataset_dim)
    with open(f'{output_dir}/.ALLCools', 'w') as f:
        yaml.dump(config, f)
    return

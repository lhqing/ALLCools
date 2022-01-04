import yaml


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

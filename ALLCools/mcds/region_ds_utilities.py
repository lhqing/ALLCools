import yaml


def update_region_ds_config(
    output_dir, new_dataset_dim=None, change_region_dim=None, config=None
):
    # update RegionDS default dimension
    try:
        with open(f"{output_dir}/.ALLCools", "r") as f:
            _config = yaml.load(f, yaml.SafeLoader)
    except FileNotFoundError:
        _config = {"region_dim": None, "ds_region_dim": {}}

    if config is not None:
        _config.update(config)

    if change_region_dim is not None:
        _config["region_dim"] = change_region_dim

    if new_dataset_dim is not None:
        _config["ds_region_dim"].update(new_dataset_dim)

    with open(f"{output_dir}/.ALLCools", "w") as f:
        yaml.dump(_config, f)
    return

import pathlib

import yaml

import ALLCools


def load_yaml(yaml_path):
    with open(yaml_path) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    if config is None:
        config = {}
    return config


# load default yaml config
DEFAULT_CONFIG_PATH = pathlib.Path(ALLCools.__path__[0]) / "dataset/dataset.yaml"
DEFAULT_CONFIG = load_yaml(DEFAULT_CONFIG_PATH)


class DatasetDirMixIn:
    def __init__(self, dataset_dir, config=None):
        self.dataset_dir = pathlib.Path(dataset_dir).absolute().resolve()

        # deal with config
        self.config = DEFAULT_CONFIG
        self._get_user_config(config)

        # deal with cell metadata
        self.metadata_dir = self.dataset_dir / self.config["meta_dir"]
        self.metadata_path = self.dataset_dir / self.config["cell_metadata"]

        # deal with mcds paths
        self.mcds_dir = self.dataset_dir / self.config["mcds_dir"]
        self.mcds_paths = []
        self._get_mcds_paths()
        return

    def _get_user_config(self, config):
        if config is not None:
            if pathlib.Path(str(config)).exists():
                config = load_yaml(config)
            elif self.dataset_dir.joinpath(config).exists():
                config = load_yaml(self.dataset_dir.joinpath(config))
        elif self.dataset_dir.joinpath("dataset.yaml").exists():
            config = load_yaml(self.dataset_dir.joinpath("dataset.yaml"))
        else:
            # no user config found, use default config
            config = {}
        self.config.update(config)
        return

    def _get_mcds_paths(self):
        all_flat_paths = []

        mcds_paths = self.config["mcds_paths"]
        if isinstance(mcds_paths, str):
            mcds_paths = [mcds_paths]

        for mcds_path in mcds_paths:
            if "*" in mcds_path:
                all_flat_paths += self.dataset_dir.glob(mcds_path)
            else:
                all_flat_paths.append(self.dataset_dir.joinpath(mcds_path))

        self.mcds_paths = all_flat_paths
        return

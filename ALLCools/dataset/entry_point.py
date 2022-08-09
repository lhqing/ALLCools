from .directory_structure import DatasetDirMixIn


class ALLCoolsDataset(DatasetDirMixIn):
    def __init__(self, dataset_dir):
        super().__init__(dataset_dir)
        return

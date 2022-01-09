import pandas as pd
import xarray as xr
import numpy as np
import pathlib
import yaml
from scipy.sparse import coo_matrix


def _make_hypo_hyper_matrix(series, dmr_values):
    samples = dmr_values.columns
    sample_int = {sample: i for i, sample in enumerate(samples)}
    rows = []
    cols = []
    datas = []
    for i, (_, samples) in enumerate(series.dropna().iteritems()):
        for sample in samples.split(","):
            rows.append(i)
            cols.append(sample_int[sample])
            datas.append(1)
    matrix = coo_matrix(
        (datas, (rows, cols)), shape=dmr_values.shape, dtype=np.int16
    ).toarray()
    matrix = pd.DataFrame(matrix, index=dmr_values.index, columns=dmr_values.columns)
    return matrix


def methylpy_to_region_ds(dmr_path, output_dir):
    pathlib.Path(output_dir).mkdir(parents=True)
    with open(f"{output_dir}/.ALLCools", "w") as f:
        config = {"region_dim": "dmr", "ds_region_dim": {"dmr": "dmr"}}
        yaml.dump(config, f)

    methylpy_dmr = pd.read_csv(dmr_path, sep="\t")

    # process index
    dmr_name = "_".join(dmr_path.split("/")[-1].split("_")[:-3])
    methylpy_dmr.index = methylpy_dmr.index.map(lambda i: f"{dmr_name}-{i}")
    methylpy_dmr.index.name = "dmr"
    methylpy_dmr.columns.name = "sample"

    dmr_infos = methylpy_dmr.iloc[:, :6]
    dmr_values = methylpy_dmr.iloc[:, 6:]
    dmr_values.columns = dmr_values.columns.map(lambda i: "_".join(i.split("_")[2:]))

    hyper_matrix = _make_hypo_hyper_matrix(
        dmr_infos["hypermethylated_samples"], dmr_values
    )
    hypo_matrix = _make_hypo_hyper_matrix(
        dmr_infos["hypomethylated_samples"], dmr_values
    )
    dmr_state = hyper_matrix - hypo_matrix
    dmr_ds = xr.Dataset(
        {"dmr_state": dmr_state, "dmr_da_frac": xr.DataArray(dmr_values)}
    )

    dmr_ds.coords.update(
        {
            "dmr_chrom": dmr_infos.iloc[:, 0],
            "dmr_start": dmr_infos.iloc[:, 1],
            "dmr_end": dmr_infos.iloc[:, 2] + 2,
            "dmr_ndms": dmr_infos.iloc[:, 3],
        }
    )

    dmr_ds.coords["dmr_chrom"] = dmr_ds.coords["dmr_chrom"].astype("str")
    dmr_ds.coords["sample"] = dmr_ds.coords["sample"].astype("str")
    dmr_ds.coords["dmr"] = dmr_ds.coords["dmr"].astype("str")

    dmr_ds["dmr_state"] = dmr_ds["dmr_state"].transpose("sample", "dmr")
    dmr_ds["dmr_da_frac"] = dmr_ds["dmr_da_frac"].transpose("sample", "dmr")

    dmr_ds.to_zarr(f"{output_dir}/dmr")
    return

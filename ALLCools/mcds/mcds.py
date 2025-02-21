import glob
import pathlib
import re
import warnings
from typing import Union

import anndata
import dask
import numpy as np
import pandas as pd
import scipy.sparse as ss
import xarray as xr
import yaml
from pybedtools import BedTool

from ..plot.qc_plots import cutoff_vs_cell_remain, plot_dispersion
from .utilities import (
    calculate_posterior_mc_frac,
    determine_engine,
    highly_variable_methylation_feature,
)


def _make_obs_df_var_df(use_data, obs_dim, var_dim):
    obs_df = pd.DataFrame([], index=use_data.get_index(obs_dim).astype(str))
    var_df = pd.DataFrame([], index=use_data.get_index(var_dim).astype(str))
    coord_prefix = re.compile(f"({obs_dim}|{var_dim})_")
    for k, v in use_data.coords.items():
        if k in [obs_dim, var_dim]:
            continue
        try:
            # v.dims should be size 1
            if v.dims[0] == obs_dim:
                series = v.to_pandas()
                # adata.obs_name is str type
                series.index = series.index.astype(str)
                obs_df[coord_prefix.sub("", k)] = series
            elif v.dims[0] == var_dim:
                series = v.to_pandas()
                # adata.var_name is str type
                series.index = series.index.astype(str)
                var_df[coord_prefix.sub("", k)] = series
            else:
                pass
        except IndexError:
            # v.dims is 0, just ignore
            pass
    return obs_df, var_df


class MCDS(xr.Dataset):
    """The MCDS Class."""

    __slots__ = ()

    def __init__(self, dataset, coords=None, attrs=None, obs_dim=None, var_dim=None):
        if isinstance(dataset, xr.Dataset):
            data_vars = dataset.data_vars
            coords = dataset.coords if coords is None else coords
            attrs = dataset.attrs if attrs is None else attrs
        else:
            data_vars = dataset
        super().__init__(data_vars=data_vars, coords=coords, attrs=attrs)
        self.obs_dim = obs_dim
        self.var_dim = var_dim

        if self.obs_dim is not None:
            # validate obs_dim is unique
            if self.obs_names.duplicated().sum() != 0:
                print("Warning: obs_names are not unique.")
        return

    @property
    def var_dim(self):
        """Name of the feature dimension."""
        if self.attrs["var_dim"] == "null":
            return None
        else:
            return self.attrs["var_dim"]

    @var_dim.setter
    def var_dim(self, var_dim):
        """Set var_dim."""
        if var_dim is not None:
            if var_dim not in self.dims:
                raise KeyError(f"{var_dim} does not occur in dimension names: {list(self.dims.keys())}")
        else:
            var_dim = "null"
        self.attrs["var_dim"] = var_dim
        return

    @property
    def obs_dim(self):
        """Name of the observation dimension."""
        if self.attrs["obs_dim"] == "null":
            return None
        else:
            return self.attrs["obs_dim"]

    @obs_dim.setter
    def obs_dim(self, obs_dim):
        """Set obs_dim."""
        if obs_dim is not None:
            if obs_dim not in self.dims:
                raise KeyError(f"{obs_dim} does not occur in dimension names: {list(self.dims.keys())}")
        else:
            obs_dim = "null"
        self.attrs["obs_dim"] = obs_dim
        return

    @property
    def obs_names(self):
        """Get obs_names."""
        if self.obs_dim is None:
            return None
        else:
            obs_index = self.get_index(self.obs_dim)
            return obs_index

    @property
    def var_names(self):
        """Get var_names."""
        if self.var_dim is None:
            return None
        else:
            var_index = self.get_index(self.var_dim)
            return var_index

    def _verify_dim(self, dim, mode):
        if dim is None:
            if mode == "var":
                if self.var_dim is None:
                    raise ValueError("MCDS does not have var_dim set, please provide the name of var_dim")
                else:
                    return self.var_dim
            else:
                if self.obs_dim is None:
                    raise ValueError("MCDS does not have obs_dim set, please provide the name of obs_dim")
                else:
                    return self.obs_dim
        else:
            if dim not in self.dims:
                raise KeyError(f"{dim} does not occur in dimension names: {list(self.dims.keys())}")
            return dim

    @classmethod
    def get_var_dims(cls, mcds_paths):
        """Get var_dim from MCDS files."""
        _var_dims = set()
        if isinstance(mcds_paths, (str, pathlib.Path)):
            mcds_paths = [mcds_paths]
        for mcds_path in mcds_paths:
            config_path = pathlib.Path(f"{mcds_path}/.ALLCools")
            if config_path.exists():
                with open(config_path) as f:
                    config = yaml.safe_load(f)
                    for dim in config["ds_region_dim"].values():
                        _var_dims.add(dim)
        return list(_var_dims)

    @classmethod
    def open(
        cls,
        mcds_paths,
        obs_dim="cell",
        use_obs=None,
        var_dim=None,
        chunks="auto",
        split_large_chunks=False,
        obj_to_str=True,
        engine=None,
        coords="minimal",
        compat="override",
        **kwargs,
    ):
        """
        Take one or multiple MCDS file paths and create single MCDS concatenated on obs_dim.

        Parameters
        ----------
        mcds_paths
            Single MCDS path or MCDS path pattern with wildcard or MCDS path list
        obs_dim
            Dimension name of observations, default is 'cell'
        use_obs
            Subset the MCDS by a list of observation IDs.
        var_dim
            Which var_dim dataset to use, needed when MCDS has multiple var_dim stored in the same directory
        chunks
            if not None, xarray will use chunks to load data as dask.array. The "auto" means xarray will
            determine chunks automatically. For more options, read the `xarray.open_dataset` `chunks` parameter
            documentation. If None, xarray will not use dask, which is not desired in most cases.
        split_large_chunks
            Whether split large chunks in dask config array.slicing.split_large_chunks
        obj_to_str
            Whether turn object coordinates into string data type
        engine
            xarray engine used to store MCDS, if multiple MCDS provided, the engine need to be the same
        coords
            the coords parameter of :py:func:`xarray.open_mfdataset` function,
            default is "minimal", means only coordinates in which the dimension already appears are included.
        compat
            the compat parameter of :py:func:`xarray.open_mfdataset` function,
            default is "override", means skip comparing variables with the same name and pick variable from first MCDS.
        kwargs
            Additional arguments passed on to :py:func:`xarray.open_dataset`
            or :py:func:`xarray.open_mfdataset` function.

        Returns
        -------
        MCDS
        """
        # parse wildcard, if any
        if isinstance(mcds_paths, (str, pathlib.Path)):
            mcds_paths = [mcds_paths]
        _flat_mcds_paths = []
        for path in mcds_paths:
            if isinstance(path, str):
                if "*" in path:
                    _flat_mcds_paths += list(glob.glob(path))
                else:
                    _flat_mcds_paths.append(path)
            else:
                _flat_mcds_paths.append(str(path))

        # if var_dim is a list, open each separately and return merged MCDS:
        if isinstance(var_dim, list):
            ds_list = []
            for _var_dim in var_dim:
                ds = cls.open(
                    _flat_mcds_paths,
                    obs_dim=obs_dim,
                    use_obs=use_obs,
                    var_dim=_var_dim,
                    chunks=chunks,
                    split_large_chunks=split_large_chunks,
                    engine=engine,
                    **kwargs,
                )
                ds_list.append(ds)
            ds = xr.merge(ds_list)
            return cls(ds, obs_dim=obs_dim, var_dim=None).squeeze()

        # determine dataset var_dim
        _var_dims = set()
        _final_paths = []
        has_dataset = False
        for path in _flat_mcds_paths:
            config_path = pathlib.Path(f"{path}/.ALLCools")
            if config_path.exists():
                has_dataset = True
                with open(config_path) as f:
                    config = yaml.safe_load(f)
                    if var_dim is None:
                        if config["region_dim"] is None:
                            # try to infer if there is only one region
                            if len(config["ds_region_dim"]) == 1:
                                _ds_var_dim = list(config["ds_region_dim"].values())[0]
                                _var_dims.add(_ds_var_dim)
                            else:
                                raise ValueError(
                                    f"MCDS paths containing multiple datasets "
                                    f'{list(config["ds_region_dim"].keys())}, '
                                    'please specify one or multiple names with the "var_dim" parameter.'
                                )
                        else:
                            _ds_var_dim = config["region_dim"]
                            _var_dims.add(config["region_dim"])
                        _final_paths.append(f"{path}/{_ds_var_dim}")
                    else:
                        if var_dim not in config["ds_region_dim"].keys():
                            raise KeyError(f"{path} do not have {var_dim}")
                        _var_dims.add(var_dim)
                        _final_paths.append(f"{path}/{var_dim}")
            else:
                # single netcdf or zarr file
                _final_paths.append(path)

        if has_dataset:
            if len(_var_dims) != 1:
                raise ValueError("Some MCDS dataset has multiple var_dim, please specify var_dim parameter.")
            else:
                var_dim = _var_dims.pop()
        for path in _final_paths:
            if not pathlib.Path(path).exists():
                raise FileNotFoundError(f"{path} does not exist.")

        # determine engine
        if engine is None:
            engine = determine_engine(_final_paths)

        if len(_final_paths) == 1:
            ds = xr.open_dataset(
                _final_paths[0],
                engine=engine,
                **kwargs,
                # chunks=chunks  # do not apply chunks parameter here, apply in the end
            )
        else:
            with dask.config.set(**{"array.slicing.split_large_chunks": split_large_chunks}):
                kwargs["coords"] = coords
                kwargs["compat"] = compat
                ds = xr.open_mfdataset(
                    _final_paths,
                    parallel=False,
                    combine="nested",
                    concat_dim=obs_dim,
                    # chunks=chunks,  # do not apply chunks parameter here, apply in the end
                    engine=engine,
                    **kwargs,
                )

        if use_obs is not None:
            with dask.config.set(**{"array.slicing.split_large_chunks": split_large_chunks}):
                use_obs_bool = ds.get_index(obs_dim).isin(use_obs)
                ds = ds.sel({obs_dim: use_obs_bool})

        # change object dtype to fix-sized string dtype
        if obj_to_str:
            for k in ds.coords.keys():
                if ds.coords[k].dtype == "O":
                    # must load then convert dtype,
                    # otherwise the size somehow not converted correctly
                    ds.coords[k] = ds.coords[k].load().astype(str)

        if chunks is not None:
            if chunks == "auto":
                n_obs = ds.get_index(obs_dim).size
                data_chunks = {}
                for name, da in ds.data_vars.items():
                    if "chunks" in da.encoding:
                        # if encoding exist, use the chunks in encoding
                        data_chunks[name] = da.encoding["chunks"]
                    else:
                        # else, use the default chunk setting
                        chunks = {
                            obs_dim: max(min(1000, max(n_obs // 10, 1)), 10),
                            # by default, no chunk on features
                            "mc_type": 1,
                            "count_type": 1,
                        }
                        data_chunks[name] = tuple(chunks[k] if k in chunks else da.get_index(k).size for k in da.dims)
                for name, da in ds.data_vars.items():
                    if da.chunks is None:
                        ds[name] = da.chunk(data_chunks[name])
            else:
                # user provide specific chunks
                if isinstance(chunks, dict):
                    _chunks = {dim: chunks[dim] if dim in chunks else ds.get_index(dim).size for dim in ds.dims}
                    ds = ds.chunk(chunks=_chunks)
                else:
                    raise ValueError(f'chunks must be a dict or "auto", got {chunks}')

        return cls(ds, obs_dim=obs_dim, var_dim=var_dim).squeeze()

    def add_mc_frac(
        self,
        var_dim=None,
        da=None,
        normalize_per_cell=True,
        clip_norm_value=10,
        da_suffix="frac",
    ):
        """
        Add posterior mC rate data array for certain feature type (var_dim).

        Parameters
        ----------
        var_dim
            Name of the feature type
        da
            if None, will use f'{var_dim}_da'
        normalize_per_cell
            if True, will normalize the mC rate data array per cell
        clip_norm_value
            reset larger values in the normalized mC rate data array to this
        da_suffix
            name suffix appended to the calculated mC rate data array
        """
        var_dim = self._verify_dim(var_dim, mode="var")

        if da is None:
            if var_dim is None:
                raise ValueError("Both da and var_dim are not specified")
            da = f"{var_dim}_da"

        if da not in self.data_vars:
            raise KeyError(f"{da} is not in this dataset")
        if var_dim not in self[da].dims:
            raise KeyError(f"{var_dim} is not a dimension of {da}")

        frac = self._calculate_frac(
            var_dim=var_dim, da=da, normalize_per_cell=normalize_per_cell, clip_norm_value=clip_norm_value
        )
        self[da + "_" + da_suffix] = frac
        return

    def _calculate_frac(self, var_dim, da, normalize_per_cell, clip_norm_value):
        """Calculate mC frac data array for certain feature type (var_dim)."""
        var_dim = self._verify_dim(dim=var_dim, mode="var")

        da_mc = self[da].sel(count_type="mc")
        da_cov = self[da].sel(count_type="cov")
        frac = calculate_posterior_mc_frac(
            mc_da=da_mc,
            cov_da=da_cov,
            var_dim=var_dim,
            normalize_per_cell=normalize_per_cell,
            clip_norm_value=clip_norm_value,
        )
        return frac

    def add_m_value(self, var_dim=None, da=None, alpha=0.01, da_suffix="mvalue"):
        """
        Add m value data array for certain feature type (var_dim).

        M-Value is a transformation of the posterior mC fraction data array to a log ratio scale.
        M = np.log2((frac + alpha) / (1 - frac + alpha)).

        No normalization is applied to M-Value.

        Parameters
        ----------
        var_dim :
            Name of the feature type
        da :
            DataArray name. if None, will use f'{var_dim}_da'
        alpha :
            alpha value for the transformation regularization
        da_suffix :
            name suffix appended to the calculated mC rate data array
        """
        var_dim = self._verify_dim(var_dim, mode="var")

        if da is None:
            if var_dim is None:
                raise ValueError("Both da and var_dim are not specified")
            da = f"{var_dim}_da"
        if da not in self.data_vars:
            raise KeyError(f"{da} is not in this dataset, please specify da or var_dim")
        if var_dim not in self[da].dims:
            raise KeyError(f"{var_dim} is not a dimension of {da}")

        frac = self._calculate_frac(var_dim=var_dim, da=da, normalize_per_cell=False, clip_norm_value=None)

        m_value = np.log2((frac + alpha) / (1 - frac + alpha))

        self[da + "_" + da_suffix] = m_value
        return

    def add_mc_rate(self, *args, **kwargs):
        """Add mC fraction data array (Deprecated)."""
        warnings.warn(
            'MCDS.add_mc_rate is renamed to MCDS.add_mc_frac, the default suffix also changed from "rate" to "frac"',
            DeprecationWarning,
        )
        self.add_mc_frac(*args, **kwargs)

    def add_feature_cov_mean(self, obs_dim=None, var_dim=None, plot=True, da_name=None):
        """
        Add feature cov mean across obs_dim.

        Parameters
        ----------
        var_dim
            Name of var dimension
        obs_dim
            Name of obs dimension
        plot
            If true, plot the distribution of feature cov mean
        da_name
            Name of the calculated data array, if None, will use f'{var_dim}_da'
        """
        obs_dim = self._verify_dim(obs_dim, mode="obs")
        var_dim = self._verify_dim(var_dim, mode="var")

        if da_name is None:
            da_name = f"{var_dim}_da"

        cov_mean_key = f"{var_dim}_cov_mean"
        if cov_mean_key in self.coords:
            print(
                f"Using existing {cov_mean_key} in MCDS.coords, if this is not desired, delete the "
                f'existing results by `del mcds.coords["{cov_mean_key}"]` and recalculate.'
            )
            feature_cov_mean = self.coords[cov_mean_key].to_pandas()
        else:
            _da = self[da_name]
            if "mc_type" in _da.dims:
                feature_cov_mean = _da.sel(count_type="cov").sum(dim="mc_type").mean(dim=obs_dim).squeeze().to_pandas()
            else:
                feature_cov_mean = _da.sel(count_type="cov").mean(dim=obs_dim).squeeze().to_pandas()
            self.coords[cov_mean_key] = feature_cov_mean
            print(f"Feature {var_dim} mean cov across cells added in MCDS.coords['{cov_mean_key}'].")
        if plot:
            cutoff_vs_cell_remain(feature_cov_mean, name=f"{var_dim}_cov_mean")
        return

    def add_cell_metadata(self, metadata, obs_dim=None):
        """Add cell metadata table to the MCDS."""
        obs_dim = self._verify_dim(obs_dim, mode="obs")
        metadata.index.name = obs_dim
        mcds_index = self.get_index(obs_dim)
        for name, col in metadata.reindex(mcds_index).items():
            self.coords[f"{obs_dim}_{name}"] = col
        return

    def filter_feature_by_cov_mean(self, var_dim=None, min_cov=0, max_cov=999999):
        """
        Filter MCDS by feature cov mean. add_feature_cov_mean() must be called before this function.

        Parameters
        ----------
        var_dim
            Name of var dimension
        min_cov
            Minimum cov cutoff
        max_cov
            Maximum cov cutoff

        Returns
        -------
        MCDS
        """
        var_dim = self._verify_dim(var_dim, mode="var")

        try:
            feature_cov_mean = self.coords[f"{var_dim}_cov_mean"].to_pandas()
        except KeyError:
            raise KeyError(f"{var_dim}_cov_mean not found in the coords, run add_feature_cov_mean() first.")

        judge: pd.Series = (feature_cov_mean > min_cov) & (feature_cov_mean < max_cov)
        with dask.config.set(**{"array.slicing.split_large_chunks": False}):
            # ignore dask warning
            mcds = self.sel({var_dim: judge[judge].index})
        selected_feature = judge.sum()
        ratio = selected_feature / judge.size * 100
        print(f"Before cov mean filter: {judge.size} {var_dim}")
        print(f" After cov mean filter: {selected_feature} {var_dim} {ratio:.1f}%")
        return mcds

    def get_feature_bed(self, var_dim=None):
        """
        Get a bed format data frame of the var_dim.

        Parameters
        ----------
        var_dim
            Name of var_dim

        Returns
        -------
        pd.DataFrame
        """
        var_dim = self._verify_dim(var_dim, mode="var")
        if f"{var_dim}_bin_start" in self:
            start_name = f"{var_dim}_bin_start"
        elif f"{var_dim}_start" in self:
            start_name = f"{var_dim}_start"
        else:
            raise KeyError(f"Neither {var_dim}_bin_start nor {var_dim}_start exists in MCDS")
        if f"{var_dim}_bin_end" in self:
            end_name = f"{var_dim}_bin_end"
        elif f"{var_dim}_end" in self:
            end_name = f"{var_dim}_end"
        else:
            raise KeyError(f"Neither {var_dim}_bin_end nor {var_dim}_end exists in MCDS")

        bed_df = pd.DataFrame(
            [
                self.get(f"{var_dim}_chrom").to_pandas(),
                self.get(start_name).to_pandas(),
                self.get(end_name).to_pandas(),
            ],
            index=["chrom", "start", "end"],
            columns=self.get_index(var_dim),
        ).T
        return bed_df

    def remove_black_list_region(self, black_list_path, var_dim=None, f=0.2):
        """
        Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path.

        Parameters
        ----------
        var_dim
            Name of var_dim
        black_list_path
            Path to the black list bed file
        f
            Fraction of overlap when calling bedtools intersect

        Returns
        -------
        MCDS
        """
        var_dim = self._verify_dim(var_dim, mode="var")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            feature_bed_df = self.get_feature_bed(var_dim=var_dim)
            feature_bed = BedTool.from_dataframe(feature_bed_df)
            black_list_bed = BedTool(black_list_path)
            black_feature = feature_bed.intersect(black_list_bed, f=f, wa=True)
            try:
                black_feature_index = black_feature.to_dataframe().set_index(["chrom", "start", "end"]).index
            except pd.errors.EmptyDataError:
                print("No feature overlapping the black list bed file.")
                return self
        black_feature_id = pd.Index(
            feature_bed_df.reset_index().set_index(["chrom", "start", "end"]).loc[black_feature_index][var_dim]
        )

        print(
            f"{black_feature_id.size} {var_dim} features removed due to overlapping"
            f" (bedtools intersect -f {f}) with black list regions."
        )
        with dask.config.set(**{"array.slicing.split_large_chunks": False}):
            mcds = self.sel({var_dim: ~self.get_index(var_dim).isin(black_feature_id)})
        return mcds

    def remove_chromosome(self, exclude_chromosome=None, include_chromosome=None, var_dim=None):
        """
        Remove regions in specific chromosome.

        Parameters
        ----------
        var_dim
            Name of var_dim
        exclude_chromosome
            if provided, only these chromosomes will be removed
        include_chromosome
            if provided, only these chromosomes will be kept

        Returns
        -------
        MCDS (xr.Dataset)
        """
        var_dim = self._verify_dim(var_dim, mode="var")
        if exclude_chromosome is not None:
            judge = self.get(f"{var_dim}_chrom").isin(exclude_chromosome)
            print(f"{int(judge.sum())} {var_dim} features in {exclude_chromosome} removed.")
            with dask.config.set(**{"array.slicing.split_large_chunks": False}):
                mcds = self.sel({var_dim: ~judge})
        if include_chromosome is not None:
            judge = self.get(f"{var_dim}_chrom").isin(include_chromosome)
            print(f"{int(judge.sum())} {var_dim} features in {include_chromosome} kept.")
            with dask.config.set(**{"array.slicing.split_large_chunks": False}):
                mcds = self.sel({var_dim: judge})
        return mcds

    def _get_da_name(self, var_dim, da_suffix):
        name1 = f"{var_dim}_{da_suffix}"
        if name1 in self:
            return name1
        name2 = f"{var_dim}_da_{da_suffix}"
        if name2 in self:
            return name2
        raise KeyError(f"{name1} or {name2} not found in MCDS.")

    def calculate_hvf_svr(
        self, mc_type=None, var_dim=None, obs_dim=None, n_top_feature=5000, da_name=None, da_suffix="frac", plot=True
    ):
        """Calculate the highly variable features (hvf) with the Support Vector Regression model."""
        from sklearn.svm import SVR

        obs_dim = self._verify_dim(obs_dim, mode="obs")
        var_dim = self._verify_dim(var_dim, mode="var")

        if da_name is None:
            da_name = f"{var_dim}_da"
            frac_da_name = self._get_da_name(var_dim, da_suffix)
        else:
            frac_da_name = f"{da_name}_{da_suffix}"

        frac_da = self[frac_da_name]
        count_da = self[da_name]
        if "mc_type" in self.dims:
            if mc_type is None:
                mc_types = self.get_index("mc_type")
                if mc_types.size > 1:
                    raise ValueError("Data array has multiple mc_types, please specify mc_type parameter.")
                else:
                    mc_type = mc_types[0]
            feature_mc_frac_mean = frac_da.sel(mc_type=mc_type).mean(dim=obs_dim).to_pandas()
            feature_std = frac_da.sel(mc_type=mc_type).std(dim=obs_dim).to_pandas()
            feature_cov_mean = count_da.sel(mc_type=mc_type, count_type="cov").mean(dim=obs_dim).to_pandas()
        else:
            feature_mc_frac_mean = frac_da.mean(dim=obs_dim).to_pandas()
            feature_std = frac_da.std(dim=obs_dim).to_pandas()
            feature_cov_mean = count_da.sel(count_type="cov").mean(dim=obs_dim).to_pandas()

        # remove bad features
        judge = (feature_mc_frac_mean > 0) & (feature_std > 0) & (feature_cov_mean > 0)
        if n_top_feature >= judge.size:
            n_top_feature = judge.size
            print("n_top_feature is than total number of features, will use all features")
        feature_mc_frac_mean = feature_mc_frac_mean[judge]
        feature_var = feature_std[judge] ** 2  # to be consistent with previous bin-based method, use var here
        feature_cov_mean = feature_cov_mean[judge]

        # prepare data
        dispersion = feature_var / feature_mc_frac_mean
        log2_disp = np.log2(dispersion)
        log2_mc_frac_mean = np.log2(feature_mc_frac_mean)
        log2_cov_mean = np.log2(feature_cov_mean)
        x = np.vstack((log2_mc_frac_mean, log2_cov_mean)).T

        # non-linear regression predicting dispersion using mc_frac_mean and cov_mean.
        svr_gamma = 1000 / judge.sum()
        print(
            f"Fitting SVR with gamma {svr_gamma:.4f}, "
            f"predicting feature dispersion using mc_frac_mean and cov_mean."
        )
        clf = SVR(gamma=svr_gamma)
        clf.fit(x, log2_disp)
        # Score is the relative position with respect of the fitted curve
        score = log2_disp - clf.predict(x)
        selected_feature_index = score.sort_values()[-n_top_feature:].index
        # make results table
        selected_feature_index = score.sort_values()[-n_top_feature:].index
        hvf_df = pd.DataFrame(
            {
                "mean": feature_mc_frac_mean.reindex(judge.index).fillna(0),
                "dispersion": dispersion.reindex(judge.index).fillna(0),
                "cov": feature_cov_mean.reindex(judge.index).fillna(0),
                "score": score.reindex(judge.index).fillna(-100),
            }
        )
        hvf_df["feature_select"] = hvf_df.index.isin(selected_feature_index)

        print(f"Total Feature Number:     {judge.size}")
        print(
            f"Highly Variable Feature:  {selected_feature_index.size} "
            f"({(selected_feature_index.size / judge.size * 100):.1f}%)"
        )

        if plot:
            try:
                import plotly.graph_objects as go
            except ImportError:
                print("Please install plotly to enable plotting.")
            plot = False

        if plot:
            if hvf_df.shape[0] > 5000:
                plot_data = hvf_df.sample(5000)
            else:
                plot_data = hvf_df
            fig = go.Figure(
                data=[
                    go.Scatter3d(
                        x=plot_data["mean"],
                        y=plot_data["cov"],
                        z=np.log2(plot_data["dispersion"]),
                        mode="markers",
                        hoverinfo="none",
                        marker={
                            "size": 2,
                            "color": plot_data["feature_select"].map({True: "red", False: "gray"}).tolist(),
                            "opacity": 0.8,
                        },
                    )
                ]
            )
            fig.update_layout(
                scene={
                    "xaxis_title": "mC Frac. Mean",
                    "yaxis_title": "Coverage Mean",
                    "zaxis_title": "log2(Dispersion)",
                },
                margin={"r": 0, "b": 0, "l": 0, "t": 0},
            )
            fig.show()

        for name, column in hvf_df.items():
            if mc_type is None:
                self.coords[f"{var_dim}_{name}"] = column
            else:
                self.coords[f"{var_dim}_{mc_type}_{name}"] = column
        return hvf_df

    def calculate_hvf(
        self,
        mc_type=None,
        var_dim=None,
        obs_dim=None,
        min_disp=0.5,
        max_disp=None,
        min_mean=0,
        max_mean=5,
        n_top_feature=5000,
        bin_min_features=5,
        mean_binsize=0.05,
        cov_binsize=100,
        da_name=None,
        da_suffix="frac",
        plot=True,
    ):
        """
        Calculate normalized dispersion to select highly variable features.

        Parameters
        ----------
        mc_type
            Type of mC to calculate
        var_dim
            Name of variable
        obs_dim
            Name of observation, default is cell
        min_disp
            minimum dispersion for a feature to be considered
        max_disp
            maximum dispersion for a feature to be considered
        min_mean
            minimum mean for a feature to be considered
        max_mean
            maximum mean for a feature to be considered
        n_top_feature
            Top N feature to use as highly variable feature.
            If set, all the cutoff will be ignored, HDF selected based on order of normalized dispersion.
        bin_min_features
            Minimum number of features to be considered as a separate bin,
            if bellow this number, the bin will be merged to its closest bin.
        mean_binsize
            bin size to separate features across mean
        cov_binsize
            bin size to separate features across coverage
        da_name
            Name of da to use, default is None, infer from var_dim and da_suffix
        da_suffix
            Suffix to add to the name of the dataarray
        plot
            If true, will plot mean, coverage and normalized dispersion scatter plots.

        Returns
        -------
        pd.DataFrame
        """
        obs_dim = self._verify_dim(obs_dim, mode="obs")
        var_dim = self._verify_dim(var_dim, mode="var")

        if da_name is None:
            frac_da_name = self._get_da_name(var_dim, da_suffix)
        else:
            frac_da_name = f"{da_name}_{da_suffix}"

        frac_da = self[frac_da_name]
        if "mc_type" in self.dims:
            if mc_type is None:
                mc_types = self.get_index("mc_type")
                if mc_types.size > 1:
                    raise ValueError("Data array has multiple mc_types, please specify mc_type parameter.")
                else:
                    mc_type = mc_types[0]
            matrix = frac_da.sel(mc_type=mc_type).squeeze()
        else:
            matrix = frac_da.squeeze()

        feature_mean_cov = self.coords[f"{var_dim}_cov_mean"]

        hvf_df = highly_variable_methylation_feature(
            cell_by_feature_matrix=matrix,
            feature_mean_cov=feature_mean_cov,
            obs_dim=obs_dim,
            var_dim=var_dim,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
            n_top_feature=n_top_feature,
            bin_min_features=bin_min_features,
            mean_binsize=mean_binsize,
            cov_binsize=cov_binsize,
        )

        selection = hvf_df["feature_select"]
        print(f"Total Feature Number:     {selection.size}")
        print(f"Highly Variable Feature:  {selection.sum()} ({(selection.sum() / selection.size * 100):.1f}%)")

        if plot:
            plot_dispersion(
                hvf_df,
                hue="feature_select",
                zlab="dispersion_norm",
                data_quantile=(0.01, 0.99),
                save_animate_path=None,
                fig_kws=None,
            )

        for name, column in hvf_df.items():
            if mc_type is None:
                self.coords[f"{var_dim}_{name}"] = column
            else:
                self.coords[f"{var_dim}_{mc_type}_{name}"] = column
        return hvf_df

    def get_count_adata(
        self,
        da_name,
        obs_dim=None,
        var_dim=None,
        sparse=True,
        loading_chunk=10000,
        binarize_cutoff=None,
        dtype="float32",
        use_vars=None,
        split_large_chunks=False,
    ):
        """
        Convert a cell-by-feature count dataarray to an adata object.

        Parameters
        ----------
        da_name
            Name of the dataarray
        obs_dim
            Name of observation
        var_dim
            Name of variable
        sparse
            If true, will adata.X convert to sparse matrix
        loading_chunk
            Chunk size to load data in memory
        binarize_cutoff
            If not None, will binarize the dataarray with values > binarize_cutoff as 1,
            values <= binarize_cutoff as 0
        dtype
            The final dtype of the adata.X, if None, will use the dtype of the dataarray
        use_vars
            If not None, will use the specified variables as vars
        split_large_chunks
            dask array.slicing.split_large_chunks parameter

        Returns
        -------
        anndata.AnnData
        """
        with dask.config.set(**{"array.slicing.split_large_chunks": split_large_chunks}):
            obs_dim = self._verify_dim(obs_dim, mode="obs")
            var_dim = self._verify_dim(var_dim, mode="var")
            da = self[da_name].transpose(obs_dim, var_dim)

            # load matrix by chunks
            da_obs = da.get_index(obs_dim)
            total_data = []
            for chunk_start in range(0, da_obs.size, loading_chunk):
                chunk_cells = da_obs[chunk_start : chunk_start + loading_chunk]
                chunk_da = da.sel({obs_dim: chunk_cells}).load()

                if use_vars is not None:
                    chunk_da = chunk_da.sel({var_dim: use_vars})

                if sparse:
                    chunk = ss.csr_matrix(chunk_da.values)
                else:
                    chunk = chunk_da.values

                # binarize the matrix to further reduce memory size
                if binarize_cutoff is not None:
                    chunk = (chunk > binarize_cutoff).astype(np.int8)
                total_data.append(chunk)
                print(f"Loading chunk {chunk_start}-{min(chunk_start + loading_chunk, da_obs.size)}" f"/{da_obs.size}")

        # concatenate all chunks
        if sparse:
            total_data = ss.vstack(total_data)
        else:
            total_data = np.vstack(total_data)

        obs_df, var_df = _make_obs_df_var_df(da, obs_dim, var_dim)
        if use_vars is not None:
            var_df = var_df.loc[use_vars, :].copy()

        if dtype is None:
            dtype = total_data.dtype

        adata = anndata.AnnData(X=total_data, obs=obs_df, var=var_df, dtype=dtype)

        if "chrom" in adata.var:
            chroms = adata.var["chrom"].astype("category")
            from natsort import natsorted

            adata.var["chrom"] = chroms.cat.reorder_categories(natsorted(chroms.cat.categories))
        return adata

    def get_score_adata(
        self,
        mc_type,
        quant_type,
        obs_dim=None,
        var_dim=None,
        sparse=True,
        dtype="float32",
        loading_chunk=50000,
        binarize_cutoff=None,
    ):
        """
        Convert a cell-by-feature methylation score dataarray to an adata object.

        Parameters
        ----------
        mc_type
            Name of the methylation type
        quant_type
            Name of the quantification type, can be "hypo-score" or "hyper-score"
        obs_dim
            Name of observation
        var_dim
            Name of variable
        sparse
            If true, will convert adata.X to sparse matrix
        dtype
            If not None, will use the dtype of the adata.X
        loading_chunk
            Chunk size to load data in memory
        binarize_cutoff
            If not None, will binarize the dataarray with values > binarize_cutoff as 1,
            values <= binarize_cutoff as 0

        Returns
        -------
        anndata.AnnData
        """
        QUANT_TYPES = ["hypo-score", "hyper-score"]
        if quant_type.lower().startswith("hypo"):
            quant_type = "hypo-score"
        elif quant_type.lower().startswith("hyper"):
            quant_type = "hyper-score"
        else:
            pass
        if quant_type not in QUANT_TYPES:
            raise ValueError(f"quant_type need to be in {QUANT_TYPES}, got {quant_type}.")

        obs_dim = self._verify_dim(obs_dim, mode="obs")
        var_dim = self._verify_dim(var_dim, mode="var")
        da_name = f"{var_dim}_da_{mc_type}-{quant_type}"
        adata = self.get_count_adata(
            da_name=da_name,
            obs_dim=obs_dim,
            var_dim=var_dim,
            sparse=sparse,
            dtype=dtype,
            loading_chunk=loading_chunk,
            binarize_cutoff=binarize_cutoff,
        )
        return adata

    def add_feature_selection_column(self, feature_select, col_name="VAR_DIM_feature_select", var_dim=None):
        """Manually add a feature selection column to the MCDS."""
        if var_dim is None:
            var_dim = self.var_dim
        features = self.get_index(var_dim)

        if isinstance(feature_select, list):
            feature_select = np.array(feature_select)

        if feature_select.dtype == "bool":
            feature_select = pd.Series(feature_select, index=features)
        else:
            feature_select = pd.Series(features.isin(feature_select), index=features)

        if col_name.startswith("VAR_DIM"):
            col_name = col_name.replace("VAR_DIM", var_dim)
        self.coords[col_name] = feature_select
        return

    def get_adata(
        self,
        mc_type=None,
        obs_dim=None,
        var_dim=None,
        da_name=None,
        da_suffix="frac",
        select_hvf=True,
        dtype="float32",
        split_large_chunks=False,
    ):
        """
        Get anndata from MCDS mC rate matrix.

        Parameters
        ----------
        mc_type
            mC rate type
        var_dim
            Name of variable
        da_suffix
            Suffix of mC rate matrix
        obs_dim
            Name of observation
        select_hvf
            Select HVF or not, if True, will use mcds.coords['{var_dim}_{mc_type}_feature_select'] to select HVFs
        dtype
            data type of adata.X
        split_large_chunks
            Whether split large chunks in dask config array.slicing.split_large_chunks

        Returns
        -------
        anndata.Anndata
        """
        obs_dim = self._verify_dim(obs_dim, mode="obs")
        var_dim = self._verify_dim(var_dim, mode="var")

        if "mc_type" in self.dims:
            if mc_type is None:
                mc_types = self.get_index("mc_type")
                if mc_types.size > 1:
                    raise ValueError("Data array has multiple mc_types, please specify mc_type parameter.")
                else:
                    mc_type = mc_types[0]
        else:
            mc_type = None

        if da_name is None:
            frac_da_name = self._get_da_name(var_dim, da_suffix)
        else:
            frac_da_name = f"{da_name}_{da_suffix}"

        with dask.config.set(**{"array.slicing.split_large_chunks": split_large_chunks}):
            frac_da = self[frac_da_name]
            if mc_type is None:
                if select_hvf:
                    if isinstance(select_hvf, str):
                        col_name = select_hvf
                    else:
                        col_name = f"{var_dim}_feature_select"
                    try:
                        use_features = self.coords[col_name].to_pandas().dropna().astype(bool)
                        use_features = use_features[use_features].index
                        use_data = frac_da.sel({var_dim: use_features}).squeeze()
                    except KeyError:
                        print(
                            f"select_hvf is {select_hvf}, but no highly variable feature results found, "
                            "use all features instead."
                        )
                        use_data = frac_da.squeeze()
                else:
                    use_data = frac_da.squeeze()
            else:
                if select_hvf:
                    if isinstance(select_hvf, str):
                        col_name = select_hvf
                    else:
                        col_name = f"{var_dim}_{mc_type}_feature_select"
                    try:
                        use_features = self.coords[col_name].to_pandas().dropna().astype(bool)
                        use_features = use_features[use_features].index
                        use_data = frac_da.sel({"mc_type": mc_type, var_dim: use_features}).squeeze()
                    except KeyError:
                        print(
                            f"select_hvf is {select_hvf}, but no highly variable feature results found, "
                            "use all features instead."
                        )
                        use_data = frac_da.sel({"mc_type": mc_type}).squeeze()
                else:
                    use_data = frac_da.sel({"mc_type": mc_type}).squeeze()

            obs_df, var_df = _make_obs_df_var_df(use_data, obs_dim, var_dim)

            adata = anndata.AnnData(
                X=use_data.astype(dtype).transpose(obs_dim, var_dim).values, obs=obs_df, var=var_df, dtype=dtype
            )
        return adata

    def merge_cluster(
        self, cluster_col, obs_dim=None, add_mc_frac=True, add_overall_mc=True, overall_mc_da="chrom100k_da"
    ):
        """Merge cell MCDS into cluster MCDS by sum on the obs_dim."""
        obs_dim = self._verify_dim(obs_dim, mode="obs")

        if isinstance(cluster_col, str):
            cluster_col = cluster_col
        else:
            self.coords[cluster_col.name] = cluster_col
            cluster_col = cluster_col.name

        with dask.config.set(**{"array.slicing.split_large_chunks": False}):
            cluster_mcds = self.groupby(cluster_col).sum(dim=obs_dim).load()
            cluster_mcds = MCDS(cluster_mcds)

            if add_mc_frac:
                for name, da in cluster_mcds.data_vars.items():
                    if name.endswith("_da") and ("count_type" in da.dims):
                        cluster_mcds.add_mc_frac(var_dim=name[:-3])

            if add_overall_mc:
                if "mc_type" in cluster_mcds.dims:
                    for mc_type in cluster_mcds.get_index("mc_type"):
                        overall_mc_dim = overall_mc_da[:-3]
                        mc = cluster_mcds[overall_mc_da].sel(mc_type=mc_type, count_type="mc").sum(dim=overall_mc_dim)
                        cov = cluster_mcds[overall_mc_da].sel(mc_type=mc_type, count_type="cov").sum(dim=overall_mc_dim)
                        cluster_mcds.coords[f"{cluster_col}_{mc_type}_overall"] = mc / cov
                else:
                    overall_mc_dim = overall_mc_da[:-3]
                    mc = cluster_mcds[overall_mc_da].sel(count_type="mc").sum(dim=overall_mc_dim)

                    cov = cluster_mcds[overall_mc_da].sel(count_type="cov").sum(dim=overall_mc_dim)
                    cluster_mcds.coords[f"{cluster_col}_overall"] = mc / cov
        return cluster_mcds

    def to_region_ds(self, region_dim=None):
        """Turn the MCDS into a RegionDS."""
        region_dim = self._verify_dim(region_dim, mode="var")

        from .region_ds import RegionDS

        _region_ds = RegionDS(self, region_dim=region_dim)
        return _region_ds

    def write_dataset(
        self, output_path, mode="w-", obs_dim=None, var_dims: Union[str, list] = None, use_obs=None, chunks="auto"
    ):
        """
        Write MCDS into an on-disk zarr dataset.

        Data arrays for each var_dim will be saved in separate sub-directories of output_path.
        The use_obs can be used to select and order observation accordingly.

        Parameters
        ----------
        output_path
            Path of the zarr dataset
        mode
            'w-' means write to output_path, fail if the path exists; 'w' means write to output_path,
            overwrite if the var_dim sub-directory exists
        obs_dim
            dimension name of observations
        var_dims
            dimension name, or a list of dimension names of variables
        use_obs
            Select and order observations according to this parameter when wrote to output_path.
        chunks
            Zarr chunks on disk and in memory.

        Returns
        -------
        output_path
        """
        with dask.config.set(**{"array.slicing.split_large_chunks": False}):
            self._write_dataset(
                output_path=output_path, mode=mode, obs_dim=obs_dim, var_dims=var_dims, use_obs=use_obs, chunks=chunks
            )
        return

    def _write_dataset(
        self, output_path, mode="w-", obs_dim=None, var_dims: Union[str, list] = None, use_obs=None, chunks="auto"
    ):
        # save mcds to a location, if the location is a MCDS dataset, auto save to sub_dir and modify .ALLCools
        # deal with large save issue, load by chunk and write append
        obs_dim = self._verify_dim(obs_dim, mode="obs")
        if var_dims is None:
            var_dims = [self._verify_dim(None, mode="var")]
        else:
            if isinstance(var_dims, str):
                var_dims = [var_dims]
            else:
                pass

        # determine output_path type
        output_path = pathlib.Path(output_path)
        config_path = pathlib.Path(f"{output_path}/.ALLCools")
        if config_path.exists():
            if mode == "w-":
                raise FileExistsError(
                    f'{output_path} already exist. Use mode="a" '
                    f"if you want to add additional data array to exist dataset. "
                    f'Or use mode="w" to overwrite existing dataset'
                )
            else:
                import shutil

                shutil.rmtree(output_path)

        # separate data_vars by var_dim
        data_vars_dict = {k: [] for k in var_dims}
        for da_name, da in self.data_vars.items():
            var_dim_count = 0
            for var_dim in var_dims:
                if var_dim in da.dims:
                    var_dim_count += 1
                    if var_dim_count == 1:
                        data_vars_dict[var_dim].append(da_name)
                    else:
                        raise ValueError(f"{da_name} contain more then one var_dim. " f"Dims of {da_name} is {da.dims}")

        # separate coordinates by var_dim
        # coords_dict = {k: [] for k in var_dims}
        # for coord_name, coord in self.coords.items():
        #     var_dim_count = 0
        #     for var_dim in var_dims:
        #         if var_dim in coord.dims:
        #             var_dim_count += 1
        #             if var_dim_count == 1:
        #                 coords_dict[var_dim].append(coord_name)
        #             else:
        #                 raise ValueError(f"{coord_name} contain more then one var_dim. "
        #                                  f"Dims of {coord_name} is {coord.dims}")

        obs_index = self.get_index(obs_dim)
        if use_obs is not None:
            # check use_obs first
            if (~use_obs.isin(obs_index)).sum() > 0:
                raise IndexError(
                    f"Not all obs id in use_obs parameter exist in MCDS "
                    f"{self.obs_dim if obs_dim is None else obs_dim} dim"
                )

        # determine chunks
        n_obs = obs_index.size
        if (chunks == "auto") or (chunks is None):
            chunks = {obs_dim: min(10000, max(n_obs // 10, 1000)), "mc_type": 1, "count_type": 1}
        full_chunks = {
            dim: size if ((dim not in chunks) or (chunks[dim] is None)) else chunks[dim]
            for dim, size in self.dims.items()
        }

        for var_dim, da_names in data_vars_dict.items():
            print(f"Saving {var_dim}")
            _this_output_path = f"{output_path}/{var_dim}"
            if pathlib.Path(_this_output_path).exists():
                print(f"{output_path}/{var_dim} will be overwrote by new data array")

            if len(da_names) == 0:
                print(f"No data for var_dim {var_dim}.")
            else:
                _mcds = self[da_names]
                # for coord_names in coords_dict.values():
                #     _mcds = _mcds.assign_coords(**{k: self.coords[k] for k in coord_names})

                # pop out obs_dim and var_dim from attrs, do not save these to dataset
                _mcds.attrs.pop("obs_dim")
                _mcds.attrs.pop("var_dim")

                total_obs = obs_index.size if use_obs is None else use_obs.size
                obs_chunk_size = chunks[obs_dim]
                for i, chunk_start in enumerate(range(0, total_obs, obs_chunk_size)):
                    print(f"Saving chunk {i}: {chunk_start} - " f"{min(chunk_start + obs_chunk_size, total_obs)}")
                    if use_obs is None:
                        _obs_chunk = obs_index[chunk_start : chunk_start + obs_chunk_size]
                        _chunk_ds = _mcds.sel({obs_dim: _obs_chunk}).load()
                    else:
                        _obs_chunk = use_obs[chunk_start : chunk_start + obs_chunk_size]
                        # load first, order next, to prevent unordered index causing more chunks
                        _chunk_ds = _mcds.sel({obs_dim: obs_index.isin(_obs_chunk)}).load()
                        _chunk_ds = _chunk_ds.sel({obs_dim: _obs_chunk})

                    # add chunk to encoding
                    for var in _chunk_ds.variables:
                        _chunk_ds[var].encoding["chunks"] = tuple(full_chunks[d] for d in _chunk_ds[var].dims)
                        if "preferred_chunks" in _chunk_ds[var].encoding:
                            _chunk_ds[var].encoding["preferred_chunks"] = {
                                d: full_chunks[d] for d in _chunk_ds[var].dims
                            }
                        if "compressor" in _chunk_ds[var].encoding:
                            try:
                                _chunk_ds[var].encoding["compressor"].clevel = 1
                            except AttributeError:
                                pass
                    if i == 0:
                        _chunk_ds.to_zarr(_this_output_path, mode="w")
                    else:
                        _chunk_ds.to_zarr(_this_output_path, append_dim=obs_dim)

        # write config
        from .utilities import update_dataset_config

        update_dataset_config(
            output_path,
            add_ds_region_dim={var_dim: var_dim for var_dim, da_list in data_vars_dict.items() if len(da_list) > 0},
            add_ds_sample_dim={var_dim: obs_dim for var_dim, da_list in data_vars_dict.items() if len(da_list) > 0},
        )
        return output_path

    def save_feature_chunk_data(
        self,
        da_name,
        output_zarr_path,
        da_suffix="_fc",
        var_dim=None,
        loading_chunk=1000,
        var_chunk_size=1,
        obs_chunk_size=500000,
        compress_level=1,
        dtype=None,
    ):
        """
        Save a data array to zarr dataset, which is chunked along the var_dim.

        This zarr dataset is useful when loading data from one or several specific
        features, such as making a gene plot.

        Parameters
        ----------
        da_name
            Name of data array to save.
        output_zarr_path
            Path to output zarr dataset.
        da_suffix
            Suffix to add to the name of the data array.
        var_dim
            Name of var_dim. If None, use self.var_dim.
        loading_chunk
            Number of var to load at a time.
        var_chunk_size
            the var_dim chunk size of the output zarr dataset.
        obs_chunk_size
            the obs_dim chunk size of the output zarr dataset.
        compress_level
            the compress level of the output zarr dataset.
        dtype
            the dtype of the output zarr dataset.
        """
        if var_dim is None:
            var_dim = self.var_dim
            var_names = self.var_names
        else:
            var_names = self.get_index(var_dim)
        obs_dim = self.obs_dim

        with dask.config.set(**{"array.slicing.split_large_chunks": False}):
            for chunk_start in range(0, var_names.size, loading_chunk):
                use_vars = var_names[chunk_start : chunk_start + loading_chunk]
                new_name = f"{da_name}{da_suffix}"
                chunk_data = xr.Dataset({new_name: self[da_name].sel({var_dim: use_vars}).load()})

                # safely change dtype, prevent overflow
                if dtype is not None:
                    try:
                        dtype_info = np.iinfo(dtype)
                    except ValueError:
                        dtype_info = np.finfo(dtype)
                    max_val = dtype_info.max
                    chunk_max = float(chunk_data[new_name].max())
                    if chunk_max > max_val:
                        print(
                            f"{da_name} max value {chunk_max} is larger than {dtype_info.max}, "
                            f"change dtype will set max value to {dtype_info.max}"
                        )
                        chunk_data[new_name] = xr.where(chunk_max > max_val, max_val, chunk_data[new_name])

                    min_val = dtype_info.min
                    chunk_min = float(chunk_data[new_name].min())
                    if chunk_min < min_val:
                        print(
                            f"{da_name} min value {chunk_min} is smaller than {dtype_info.min}, "
                            f"change dtype will set min value to {dtype_info.min}"
                        )
                        chunk_data[new_name] = xr.where(chunk_min < min_val, min_val, chunk_data[new_name])

                    # finally, change dtype
                    chunk_data[new_name] = chunk_data[new_name].astype(dtype)

                # specify the chunks of the output zarr dataset
                chunks = []
                preferred_chunks = {}
                for dim in chunk_data[new_name].dims:
                    if dim == obs_dim:
                        chunks.append(obs_chunk_size)
                        preferred_chunks[dim] = obs_chunk_size
                    elif dim == var_dim:
                        chunks.append(var_chunk_size)
                        preferred_chunks[dim] = var_chunk_size
                    elif dim in ["mc_type", "count_type"]:
                        chunks.append(1)
                        preferred_chunks[dim] = 1
                    else:
                        dim_size = chunk_data.get_index(dim).size
                        chunks.append(dim_size)
                        preferred_chunks[dim] = dim_size
                # change encoding to save the zarr in desired chunks
                chunk_data[new_name].encoding["chunks"] = chunks
                chunk_data[new_name].encoding["preferred_chunks"] = preferred_chunks

                try:
                    chunk_data[new_name].encoding["compressor"].clevel = compress_level
                except KeyError:
                    import zarr

                    chunk_data[new_name].encoding["compressor"] = zarr.Blosc(clevel=1)

                if chunk_start == 0:
                    chunk_data.to_zarr(output_zarr_path, mode="w")
                else:
                    chunk_data.to_zarr(output_zarr_path, append_dim=var_dim)
                print(f"Saved chunk {chunk_start}-{chunk_start + loading_chunk}/{var_names.size}")
        return

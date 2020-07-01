import anndata
import pandas as pd
import xarray as xr
from pybedtools import BedTool

from .utilities import highly_variable_methylation_feature, calculate_posterior_mc_rate
from ..plot.qc_plots import cutoff_vs_cell_remain, plot_dispersion


class MCDS(xr.Dataset):
    __slots__ = ()

    def __init__(self, dataset):
        super().__init__(data_vars=dataset.data_vars, coords=dataset.coords,
                         attrs=dataset.attrs)
        return

    @classmethod
    def open(cls, mcds_paths, obs_dim='cell', use_obs=None):
        """
        Take one or multiple MCDS file paths and create single MCDS concatenated on obs_dim

        Parameters
        ----------
        mcds_paths
            Single MCDS path or MCDS path pattern with wildcard or MCDS path list
        obs_dim
            Dimension name of obs, default is 'cell'
        use_obs
            Dimension name of obs, default is 'cell'
        Returns
        -------
        MCDS
        """
        if use_obs is not None:
            def preprocess(_ds):
                _use_obs = use_obs[use_obs.isin(_ds.get_index(obs_dim))]
                _ds = _ds.sel({obs_dim: _use_obs}).chunk(chunks={obs_dim: 1536})
                return _ds
        else:
            def preprocess(_ds):
                _ds = _ds.chunk(chunks={obs_dim: 1536})
                return _ds

        if isinstance(mcds_paths, str) and '*' not in mcds_paths:
            ds = xr.open_dataset(mcds_paths)
            ds = preprocess(ds)
        else:
            ds = xr.open_mfdataset(mcds_paths,
                                   preprocess=preprocess,
                                   parallel=False,
                                   combine='nested',
                                   concat_dim=obs_dim)
        return cls(ds)

    def add_mc_rate(self,
                    var_dim,
                    da=None,
                    normalize_per_cell=True,
                    clip_norm_value=10,
                    rate_da_suffix='rate'):
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
        rate_da_suffix
            name suffix appended to the calculated mC rate data array

        Returns
        -------

        """

        if da is None:
            da = f'{var_dim}_da'
        if da not in self.data_vars:
            raise KeyError(f'{da} is not in this dataset')
        if var_dim not in self[da].dims:
            raise KeyError(f'{var_dim} is not a dimension of {da}')

        da_mc = self[da].sel(count_type='mc')
        da_cov = self[da].sel(count_type='cov')

        rate = calculate_posterior_mc_rate(mc_da=da_mc,
                                           cov_da=da_cov,
                                           var_dim=var_dim,
                                           normalize_per_cell=normalize_per_cell,
                                           clip_norm_value=clip_norm_value)
        self[da + "_" + rate_da_suffix] = rate
        return

    def add_feature_cov_mean(self, var_dim, obs_dim='cell', plot=True):
        """
        Add feature cov mean across obs_dim

        Parameters
        ----------
        var_dim
            Name of var dimension
        obs_dim
            Name of obs dimension
        plot
            If true, plot the distribution of feature cov mean

        Returns
        -------

        """

        feature_cov_mean = self[f'{var_dim}_da'].squeeze().sel(
            count_type='cov').sum(dim='mc_type').mean(dim=obs_dim).to_pandas()
        self.coords[f'{var_dim}_cov_mean'] = feature_cov_mean

        print(f'Feature {var_dim} mean cov across cells added.')
        if plot:
            cutoff_vs_cell_remain(feature_cov_mean)
        return

    def filter_feature_by_cov_mean(self, var_dim, min_cov=0, max_cov=999999):
        """
        filter MCDS by feature cov mean. add_feature_cov_mean() must be called before this function.

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

        try:
            feature_cov_mean = self.coords[f'{var_dim}_cov_mean'].to_pandas()
        except KeyError:
            raise KeyError(f'{var_dim}_cov_mean not found in the coords, run add_feature_cov_mean() first.')

        judge: pd.Series = (feature_cov_mean > min_cov) & (feature_cov_mean < max_cov)
        mcds = self.sel({var_dim: judge[judge].index})

        selected_feature = judge.sum()
        ratio = selected_feature / judge.size * 100
        print(f'Before cov mean filter: {judge.size} {var_dim}')
        print(f' After cov mean filter: {selected_feature} {var_dim} {ratio:.1f}%')
        return mcds

    def get_feature_bed(self, var_dim):
        """
        Get a bed format data frame of the var_dim

        Parameters
        ----------
        var_dim
            Name of var_dim
        Returns
        -------
        pd.DataFrame
        """

        bed_df = pd.DataFrame([
            self.coords[f'{var_dim}_chrom'].to_pandas(),
            self.coords[f'{var_dim}_bin_start'].to_pandas(),
            self.coords[f'{var_dim}_bin_end'].to_pandas()
        ],
            index=['chrom', 'start', 'end'],
            columns=self.get_index(var_dim)).T
        return bed_df

    def remove_black_list_region(self, var_dim, black_list_path, f=0.2):
        """
        Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path

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
        MCDS (xr.Dataset)
        """
        feature_bed_df = self.get_feature_bed(var_dim=var_dim)
        feature_bed = BedTool.from_dataframe(feature_bed_df)

        black_list_bed = BedTool(black_list_path)
        black_feature = feature_bed.intersect(black_list_bed, f=f, wa=True)
        black_feature_index = black_feature.to_dataframe().set_index(
            ['chrom', 'start', 'end']).index
        black_feature_id = pd.Index(
            feature_bed_df.reset_index().set_index(['chrom', 'start', 'end']).loc[black_feature_index][var_dim])

        print(f'{black_feature_id.size} {var_dim} features removed due to overlapping'
              f' (bedtools intersect -f {f}) with black list regions.')

        mcds = self.sel({var_dim: ~self.get_index(var_dim).isin(black_feature_id)})
        return mcds

    def remove_chromosome(self, var_dim, exclude_chromosome):
        """
        Remove regions in specific chromosome

        Parameters
        ----------
        var_dim
            Name of var_dim
        exclude_chromosome
            Chromosome to remove

        Returns
        -------
        MCDS (xr.Dataset)
        """
        judge = self.coords[f'{var_dim}_chrom'].isin(exclude_chromosome)
        print(f'{int(judge.sum())} {var_dim} features in {exclude_chromosome} removed.')
        mcds = self.sel({var_dim: ~judge})
        return mcds

    def calculate_hvf(self,
                      mc_type,
                      var_dim,
                      obs_dim='cell',
                      min_disp=0.5,
                      max_disp=None,
                      min_mean=0,
                      max_mean=5,
                      n_top_feature=5000,
                      bin_min_features=5,
                      mean_binsize=0.05,
                      cov_binsize=100,
                      plot=True):
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
        plot
            If true, will plot mean, coverage and normalized dispersion scatter plots.

        Returns
        -------
        pd.DataFrame
        """
        matrix = self[f'{var_dim}_da_rate'].sel(mc_type=mc_type).squeeze()
        feature_mean_cov = self.coords[f'{var_dim}_cov_mean']

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
            cov_binsize=cov_binsize)

        selection = hvf_df['gene_subset']
        print(f'Total Feature Number:     {selection.size}')
        print(f'Highly Variable Feature:  {selection.sum()} ({(selection.sum() / selection.size * 100):.1f}%)')

        if plot:
            plot_dispersion(hvf_df,
                            hue='gene_subset',
                            zlab='dispersion_norm',
                            data_quantile=(0.01, 0.99),
                            save_animate_path=None,
                            fig_kws=None)
        return hvf_df

    def get_adata(self, mc_type, var_dim, rate_matrix_suffix='_da_rate', obs_dim='cell'):
        """
        Get anndata from MCDS mC rate matrix
        Parameters
        ----------
        mc_type
            mC rate type
        var_dim
            Name of variable
        rate_matrix_suffix
            Suffix of mC rate matrix
        obs_dim
            Name of observation

        Returns
        -------
        anndata.Anndata
        """
        use_data = self[f'{var_dim}{rate_matrix_suffix}'].sel({'mc_type': mc_type}).squeeze()

        obs_df = pd.DataFrame([], index=use_data.get_index(obs_dim))
        var_df = pd.DataFrame([], index=use_data.get_index(var_dim))
        for k, v in use_data.coords.items():
            if k in ['cell', var_dim]:
                continue
            try:
                # v.dims should be size 1
                if v.dims[0] == obs_dim:
                    obs_df[k] = v.to_pandas()
                elif v.dims[0] == var_dim:
                    var_df[k] = v.to_pandas()
                else:
                    pass
            except IndexError:
                # v.dims is 0, just ignore
                pass

        adata = anndata.AnnData(X=use_data.transpose(obs_dim, var_dim).values,
                                obs=obs_df,
                                var=var_df)
        return adata

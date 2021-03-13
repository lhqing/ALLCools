import anndata
import pandas as pd
import xarray as xr
import numpy as np
from pybedtools import BedTool
import dask
import re
import warnings
from .utilities import highly_variable_methylation_feature, calculate_posterior_mc_frac
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
            Dimension name of observations, default is 'cell'
        use_obs
            Subset the MCDS by a list of observation IDs.
        Returns
        -------
        MCDS
        """
        if use_obs is not None:
            def preprocess(_ds):
                _use_obs = use_obs & _ds.get_index(obs_dim)
                _ds = _ds.sel({obs_dim: _use_obs})
                return _ds
        else:
            preprocess = None

        if isinstance(mcds_paths, str) and '*' not in mcds_paths:
            ds = xr.open_dataset(mcds_paths)
            if preprocess is not None:
                ds = preprocess(ds)
        else:
            with dask.config.set(**{'array.slicing.split_large_chunks': False}):
                ds = xr.open_mfdataset(mcds_paths,
                                       preprocess=preprocess,
                                       parallel=False,
                                       combine='nested',
                                       concat_dim=obs_dim)
        return cls(ds).squeeze()

    def add_mc_frac(self,
                    var_dim,
                    da=None,
                    normalize_per_cell=True,
                    clip_norm_value=10,
                    da_suffix='frac'):
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
        frac = calculate_posterior_mc_frac(mc_da=da_mc,
                                           cov_da=da_cov,
                                           var_dim=var_dim,
                                           normalize_per_cell=normalize_per_cell,
                                           clip_norm_value=clip_norm_value)
        self[da + "_" + da_suffix] = frac

    def add_mc_rate(self, *args, **kwargs):
        warnings.warn(
            'MCDS.add_mc_rate is renamed to MCDS.add_mc_frac, the default suffix also changed from "rate" to "frac"',
            DeprecationWarning
        )
        self.add_mc_frac(*args, **kwargs)

    def add_feature_cov_mean(self, var_dim, obs_dim='cell', plot=True):
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
        Returns
        -------
        None
        """
        feature_cov_mean = self[f'{var_dim}_da'].squeeze().sel(
            count_type='cov').sum(dim='mc_type').mean(dim=obs_dim).to_pandas()
        self.coords[f'{var_dim}_cov_mean'] = feature_cov_mean

        print(f"Feature {var_dim} mean cov across cells added in MCDS.coords['{var_dim}_cov_mean'].")
        if plot:
            cutoff_vs_cell_remain(feature_cov_mean, name=f'{var_dim}_cov_mean')
        return

    def add_cell_metadata(self, metadata, obs_name='cell'):
        metadata.index.name = obs_name
        mcds_index = self.get_index(obs_name)
        for name, col in metadata.reindex(mcds_index).items():
            self.coords[f'{obs_name}_{name}'] = col
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
        with dask.config.set(**{'array.slicing.split_large_chunks': False}):
            # ignore dask warning
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
        MCDS
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
        with dask.config.set(**{'array.slicing.split_large_chunks': False}):
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
        with dask.config.set(**{'array.slicing.split_large_chunks': False}):
            mcds = self.sel({var_dim: ~judge})
        return mcds

    def calculate_hvf_svr(self, var_dim, mc_type, obs_dim='cell', n_top_feature=5000, da_suffix='frac', plot=True):
        from sklearn.svm import SVR
        import plotly.graph_objects as go

        feature_mc_frac_mean = self[f'{var_dim}_da_{da_suffix}'].sel(mc_type=mc_type).mean(
            dim=obs_dim).to_pandas()
        feature_std = self[f'{var_dim}_da_{da_suffix}'].sel(mc_type=mc_type).std(
            dim=obs_dim).to_pandas()
        feature_cov_mean = self[f'{var_dim}_da'].sel(
            mc_type=mc_type, count_type='cov').mean(dim=obs_dim).to_pandas()

        # remove bad features
        judge = (feature_mc_frac_mean > 0) & (feature_std > 0) & (feature_cov_mean > 0)
        if n_top_feature >= judge.size:
            raise ValueError('n_top_feature must be smaller than total number of features')
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
        print(f'Fitting SVR with gamma {svr_gamma:.4f}, '
              f'predicting feature dispersion using mc_frac_mean and cov_mean.')
        clf = SVR(gamma=svr_gamma)
        clf.fit(x, log2_disp)
        # Score is the relative position with respect of the fitted curve
        score = log2_disp - clf.predict(x)
        selected_feature_index = score.sort_values()[-n_top_feature:].index
        # make results table
        selected_feature_index = score.sort_values()[-n_top_feature:].index
        hvf_df = pd.DataFrame(
            {
                'mean': feature_mc_frac_mean.reindex(judge.index).fillna(0),
                'dispersion': dispersion.reindex(judge.index).fillna(0),
                'cov': feature_cov_mean.reindex(judge.index).fillna(0),
                'score': score.reindex(judge.index).fillna(-100)
            }
        )
        hvf_df['feature_select'] = hvf_df.index.isin(selected_feature_index)

        print(f'Total Feature Number:     {judge.size}')
        print(f'Highly Variable Feature:  {selected_feature_index.size} '
              f'({(selected_feature_index.size / judge.size * 100):.1f}%)')

        if plot:
            if hvf_df.shape[0] > 5000:
                plot_data = hvf_df.sample(5000)
            else:
                plot_data = hvf_df
            fig = go.Figure(data=[
                go.Scatter3d(
                    x=plot_data['mean'],
                    y=plot_data['cov'],
                    z=np.log2(plot_data['dispersion']),
                    mode='markers',
                    hoverinfo='none',
                    marker=dict(
                        size=2,
                        color=plot_data['feature_select'].map({
                            True: 'red',
                            False: 'gray'
                        }).tolist(),  # set color to an array/list of desired values
                        opacity=0.8), )
            ])
            fig.update_layout(scene=dict(xaxis_title='mC Frac. Mean',
                                         yaxis_title='Coverage Mean',
                                         zaxis_title='log2(Dispersion)'),
                              margin=dict(r=0, b=0, l=0, t=0))
            fig.show()

        for name, column in hvf_df.items():
            self.coords[f'{var_dim}_{mc_type}_{name}'] = column
        return hvf_df

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

        selection = hvf_df['feature_select']
        print(f'Total Feature Number:     {selection.size}')
        print(f'Highly Variable Feature:  {selection.sum()} ({(selection.sum() / selection.size * 100):.1f}%)')

        if plot:
            plot_dispersion(hvf_df,
                            hue='feature_select',
                            zlab='dispersion_norm',
                            data_quantile=(0.01, 0.99),
                            save_animate_path=None,
                            fig_kws=None)

        for name, column in hvf_df.items():
            self.coords[f'{var_dim}_{mc_type}_{name}'] = column
        return hvf_df

    def get_adata(self, mc_type, var_dim, da_suffix='frac', obs_dim='cell', select_hvf=True):
        """
        Get anndata from MCDS mC rate matrix
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

        Returns
        -------
        anndata.Anndata
        """
        if select_hvf:
            try:
                use_features = self.get_index(var_dim)[self.coords[f'{var_dim}_{mc_type}_feature_select']]
                use_data = self[f'{var_dim}_da_{da_suffix}'].sel({'mc_type': mc_type, var_dim: use_features}).squeeze()
            except KeyError:
                print('feature_select==True, but no highly variable feature results found, use all features instead.')
                use_data = self[f'{var_dim}_da_{da_suffix}'].sel({'mc_type': mc_type}).squeeze()
        else:
            use_data = self[f'{var_dim}_da_{da_suffix}'].sel({'mc_type': mc_type}).squeeze()

        obs_df = pd.DataFrame([], index=use_data.get_index(obs_dim).astype(str))
        var_df = pd.DataFrame([], index=use_data.get_index(var_dim).astype(str))
        coord_prefix = re.compile(f'({obs_dim}|{var_dim})_')
        for k, v in use_data.coords.items():
            if k in [obs_dim, var_dim]:
                continue
            try:
                # v.dims should be size 1
                if v.dims[0] == obs_dim:
                    series = v.to_pandas()
                    # adata.obs_name is str type
                    series.index = series.index.astype(str)
                    obs_df[coord_prefix.sub('', k)] = series
                elif v.dims[0] == var_dim:
                    series = v.to_pandas()
                    # adata.var_name is str type
                    series.index = series.index.astype(str)
                    var_df[coord_prefix.sub('', k)] = series
                else:
                    pass
            except IndexError:
                # v.dims is 0, just ignore
                pass

        adata = anndata.AnnData(X=use_data.transpose(obs_dim, var_dim).values,
                                obs=obs_df,
                                var=var_df)
        return adata

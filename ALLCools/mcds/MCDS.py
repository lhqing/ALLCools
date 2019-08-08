import numpy as np
import pandas as pd
import xarray as xr
from anndata import AnnData

from .utilities import calculate_posterior_mc_rate


class MCDS(xr.Dataset):
    def __init__(self, dataset):
        super().__init__(data_vars=dataset.data_vars, coords=dataset.coords,
                         attrs=dataset.attrs)
        return

    @classmethod
    def open(cls, paths, use_dataarrays=None, use_cells=None, sel_dict=None, squeeze=True, **kwargs):
        """
        Open one (using xarray.open_dataset) or multiple (using xarray.open_mfdataset) MCDS files.
        If provide "chunks" in kwargs or the paths is a str glob (contain *) or list,
        the return MCDS dataarray will be dask array.
        See how xarray and dask works for large dataset parallel computation:
        http://xarray.pydata.org/en/stable/dask.html

        Parameters
        ----------
        paths : str or sequence
            For single MCDS: single path string.
            For multiple MCDS: Either a string glob in the form "path/to/my/files/*.nc"
            or an explicit list of files to open.
        use_dataarrays
            String or List of dataarray name
        use_cells
        sel_dict
        squeeze
        kwargs
            Other key word arguments pass to xr.open_dataset or xr.open_mfdataset
        Returns
        -------
        MCDS
        """
        if isinstance(paths, str) and '*' not in paths:
            _this_ds = xr.open_dataset(paths, **kwargs)
        else:
            _this_ds = xr.open_mfdataset(paths, **kwargs)

        if use_dataarrays is not None:
            if isinstance(use_dataarrays, str):
                use_dataarrays = [use_dataarrays]
            da_dict = {da_name: da
                       for da_name, da in _this_ds.data_vars.items()
                       if da_name in use_dataarrays}
            _this_ds = xr.Dataset(da_dict)

        _sel_dict = {}
        if use_cells is not None:
            _sel_dict['cell'] = pd.Index(use_cells).drop_duplicates()
        if sel_dict is not None:
            _sel_dict.update(sel_dict)
        if len(_sel_dict) != 0:
            _this_ds = _this_ds.sel(**_sel_dict)

        if squeeze:
            _this_ds = _this_ds.squeeze()

        return cls(_this_ds)

    def filter_cell_cov(self, dim, da, mc_type,
                        min_cov=0, max_cov=None):
        """
        Filter cell by total cov for certain mc_type along certain dimension in certain dataarray.

        Parameters
        ----------
        dim
            region dimension to mean
        da
            dataarray to do calculation
        mc_type
            mc_type to mean
        min_cov
            minimum cov mean, suggest to plot distribution first.
        max_cov
            maximum cov mean, suggest ot plot distribution first.

        Returns
        -------

        """
        if dim not in self[da].dims:
            raise ValueError(f'{dim} is not a dimension of {da}')
        cell_mean = self[da] \
            .sel(count_type='cov', mc_type=mc_type) \
            .squeeze() \
            .mean(dim)
        if max_cov is None:
            max_cov = np.inf
        cell_max = cell_mean < max_cov
        cell_min = cell_mean > min_cov
        cell_mask = np.all(np.vstack([cell_max.values,
                                      cell_min.values]),
                           axis=0)
        select_index = self.get_index('cell')[cell_mask]
        return self.loc[dict(cell=select_index)]

    def filter_region_cov(self, dim, da, mc_type,
                          min_cov=None, max_cov=None,
                          white_list=None, black_list=None):
        """
        Filter cell by total cov for certain mc_type along certain dimension in certain dataarray.

        Parameters
        ----------
        dim
            region dimension to filter,
            Note when filtering region, mean is always performed on cell.
        da
            dataarray to do calculation
        mc_type
            mc_type to mean
        min_cov
            minimum cov mean, suggest to plot distribution first.
        max_cov
            maximum cov mean, suggest ot plot distribution first.
        white_list
            regions to keep no matter how
        black_list
            regions to remove no matter how

        Returns
        -------

        """
        if dim not in self[da].dims:
            raise ValueError(f'{dim} is not a dimension of {da}')
        region_mean = self[da] \
            .sel(count_type='cov', mc_type=mc_type) \
            .squeeze() \
            .mean('cell')
        if max_cov is None:
            max_cov = np.inf
        region_max = region_mean < max_cov
        region_min = region_mean > min_cov

        # deal with white black list
        black_list = [] if black_list is None else black_list
        white_list = [] if white_list is None else white_list
        black_set = set(black_list)
        white_set = set(white_list)
        black_and_white = black_set & white_set
        if len(black_and_white) != 0:
            raise ValueError(f'{len(black_and_white)} items in both black and white list. Check your input.')
        black_judge = self.get_index(dim).map(lambda i: i not in black_set).values
        white_judge = self.get_index(dim).map(lambda i: i in white_set).values
        region_mask = np.all(np.vstack([region_max.values,
                                        region_min.values,
                                        black_judge]), axis=0)
        region_mask = np.any(np.vstack([region_mask,
                                        white_judge]), axis=0)

        select_index = self.get_index(dim)[region_mask]
        filtered_ds = self.loc[{dim: select_index}]
        return MCDS(filtered_ds)

    def add_mc_rate(self, dim, da,
                    normalize_per_cell=True,
                    clip_norm_value=10,
                    rate_da_suffix='rate', inplace=True):
        if da not in self.data_vars:
            raise KeyError(f'{da} is not in this dataset')
        if dim not in self[da].dims:
            raise KeyError(f'{dim} is not a dimension of {da}')
        da_mc = self[da].sel(count_type='mc')
        da_cov = self[da].sel(count_type='cov')

        rate = calculate_posterior_mc_rate(mc_da=da_mc,
                                           cov_da=da_cov,
                                           var_dim=dim,
                                           normalize_per_cell=normalize_per_cell,
                                           clip_norm_value=clip_norm_value)
        if inplace:
            self[da + "_" + rate_da_suffix] = rate
            return
        else:
            return rate

    def add_gene_rate(self, dim='gene', da='gene_da', method='bayes',
                      normalize_per_cell=True, clip_norm_value=10,
                      rate_da_suffix='rate', inplace=True):
        if da not in self.data_vars:
            raise KeyError(f'{da} is not in this dataset')
        if dim not in self[da].dims:
            raise KeyError(f'{dim} is not a dimension of {da}')
        da_mc = self[da].sel(count_type='mc')
        da_cov = self[da].sel(count_type='cov')

        if method == 'bayes':
            rate = calculate_posterior_mc_rate(mc_da=da_mc,
                                               cov_da=da_cov,
                                               var_dim=dim,
                                               normalize_per_cell=normalize_per_cell,
                                               clip_norm_value=clip_norm_value)
        elif method == 'naive':
            # for gene, we just use normal rate
            rate = da_mc / da_cov

            if normalize_per_cell:
                cell_overall = da_mc.sum(dim='gene') / da_cov.sum(dim='gene')
                rate = rate / cell_overall
                if clip_norm_value is not None:
                    rate.values[np.where(rate.values > clip_norm_value)] = clip_norm_value
        else:
            raise ValueError('Method can only be "bayes" or "naive"')

        if inplace:
            self[da + "_" + rate_da_suffix] = rate
            return
        else:
            return rate

    def to_ann(self, da, var_dim, obs_dim='cell', **sel_kwargs):
        _data = self[da].sel(**sel_kwargs).squeeze().values

        if len(_data.shape) != 2:
            if ('mc_type' in self[da].dims) and ('mc_type' not in sel_kwargs):
                raise ValueError(f'The {da} dataarray have mc_type dimension, '
                                 f'you must provide mc_type for AnnData')
            else:
                raise ValueError(f'The {da} dataarray is not 2D after selection, '
                                 f'you must provide sel_kwargs to select additional dims, '
                                 f'AnnData can only be 2D.')

        index_dict = self[da].indexes
        obs_df = pd.DataFrame({k: pd.Series(v.tolist())
                               for k, v in index_dict.items()
                               if (v.name == obs_dim) and (k != obs_dim)},
                              index=index_dict[obs_dim])
        var_df = pd.DataFrame({k: pd.Series(v.tolist())
                               for k, v in index_dict.items()
                               if (v.name == var_dim) and (k != var_dim)},
                              index=index_dict[var_dim])

        return AnnData(X=_data.copy(),
                       obs=obs_df,
                       var=var_df)

    def add_ann_to_da(self, adata, var_dim=None, obs_dim='cell'):
        # columns from AnnData.obs and AnnData.var go to da.coords
        # obsm goes to new da with corresponding new dimension
        obs_df = adata.obs
        obs_df.index.name = obs_dim  # make sure coords added with "cell" index
        for col, data in obs_df.iteritems():
            self.coords[col] = data

        if var_dim is not None:
            var_df = adata.var
            var_df.index.name = var_dim  # make sure coords added with var index
            for col, data in var_df.iteritems():
                self.coords[col] = data

        for obsm_key in adata.obsm_keys():
            coord_name = obsm_key[2:]  # remove 'X_'
            obsm_data = adata.obsm[obsm_key]
            obsm_df = pd.DataFrame(obsm_data,
                                   index=adata.obs_names,
                                   columns=[f'{coord_name}_{i}' for i in range(obsm_data.shape[1])])
            obsm_df.index.name = obs_dim
            obsm_df.columns.name = coord_name
            self[coord_name + '_coord'] = obsm_df

        for varm_key in adata.varm_keys():
            coord_name = varm_key
            varm_data = adata.varm[varm_key]
            varm_df = pd.DataFrame(varm_data,
                                   index=adata.var_names,
                                   columns=[f'{coord_name}_{i}' for i in range(varm_data.shape[1])])
            varm_df.index.name = var_dim
            varm_df.columns.name = coord_name
            self[coord_name + '_coord'] = varm_df
        return

    def add_dataframe_to_coords(self, df, index_dim=None):
        # add columns to da.coords based on index and index_name
        if df.index.name not in self.dims:
            if index_dim is None:
                raise ValueError('Dataframe index name is not in dims, index_dim must be specified.')
            else:
                df.index.name = index_dim
        for col, data in df.iteritems():
            if col == index_dim:
                print(f'Found index_dim "{index_dim}" in the col, skip to prevent mess up. '
                      'If you want to add that column, rename it.')
                continue
            self.coords[col] = data
        return

    def add_dataframe_to_da(self, df, index_dim, col_dim, da_name):
        # add columns to da.coords based on index and index_name
        df.index.name = index_dim
        df.columns.name = col_dim
        self[da_name] = df
        return

    def get_cell_tidy_data(self, tsne=True, umap=True, pca=True, pc_components=4,  # Select coordinates
                           select_genes=None, gene_mc_type='CHN', add_gene_cov=True,
                           add_gene_rna=True, rna_count_type='gene', norm_rna_by_cell=True):
        # A cell tidy dataframe need:
        # - tsne, umap and pca coordinates
        # - gene rate and gene cov
        # - other cell metadata

        cell_coord_da = self.coords['cell'].coords
        records = []
        single_values = {}
        for coord in cell_coord_da.values():
            if len(coord.shape) == 0:
                # 0-d array will raise error in coord.to_series()
                single_values[coord.name] = coord.values.tolist()  # this will not give a list, but a value...
            else:
                records.append(coord.to_series())
        data = pd.DataFrame(records)
        data = data.T

        for col_name, value in single_values.items():
            data[col_name] = value

        all_dfs = [data]
        if tsne:
            if 'tsne_coord' not in self.data_vars:
                print('tSNE coords do not exist, skip')
            else:
                coord_df = self['tsne_coord'].to_pandas()
                all_dfs.append(coord_df)
        if umap:
            if 'umap_coord' not in self.data_vars:
                print('UMAP coords do not exist, skip')
            else:
                coord_df = self['umap_coord'].to_pandas()
                all_dfs.append(coord_df)
        if pca:
            if 'pca_coord' not in self.data_vars:
                print('PCA coords do not exist, skip')
            else:
                coord_df = self['pca_coord'][:, :pc_components].to_pandas()
                all_dfs.append(coord_df)

        if select_genes is not None:
            if isinstance(select_genes, str):
                select_genes = [select_genes]
            # select_genes must be gene_id used for gene coords.
            # there is not transformation and will raise if gene id not found
            gene_index = self.get_index('gene')
            for gene in select_genes:
                unknown_genes = []
                if gene not in gene_index:
                    unknown_genes.append(gene)
                if len(unknown_genes) != 0:
                    unknown_genes = ' '.join(unknown_genes)
                    raise KeyError(f'{unknown_genes} not exist in gene coords')

            # add gene mC info
            if 'gene_da_rate' not in self.data_vars:
                print('Gene rate dataarray do not exist, skip adding gene info')
            else:
                if gene_mc_type not in self['gene_da_rate'].coords['mc_type']:
                    raise ValueError(f'Gene rate dataarray do not contain mC type {gene_mc_type}')

                gene_df = self['gene_da_rate'].sel(mc_type=gene_mc_type, gene=select_genes).to_pandas()
                all_dfs.append(gene_df)
                if add_gene_cov:
                    gene_cov_df = self['gene_da'].sel(mc_type=gene_mc_type,
                                                      count_type='cov',
                                                      gene=select_genes).to_pandas()
                    gene_cov_df.columns = gene_cov_df.columns.map(lambda i: i + '_cov')
                    all_dfs.append(gene_cov_df)

            # add gene RNA info, only when rna_da exist (e.g. snmCT-seq)
            if add_gene_rna and 'rna_da' in self.data_vars:
                rna_df = self['rna_da'].sel(rna_count_type=rna_count_type,
                                            gene=select_genes).to_pandas()
                if norm_rna_by_cell:
                    # TODO optional log1p
                    # CPM, but not RPKM or TPM
                    rna_total = self['rna_da'].sel(rna_count_type=rna_count_type).sum(dim='gene')
                    rna_df = rna_df.divide(rna_total, axis=0) * 1000000

                rna_df.columns = rna_df.columns.map(lambda i: i + '_rna')
                all_dfs.append(rna_df)

        total_df = pd.concat(all_dfs, axis=1, sort=True)
        return total_df

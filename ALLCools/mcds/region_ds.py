import pandas as pd
import pybedtools
import xarray as xr
import numpy as np
import dask
import warnings
import joblib
import subprocess
import pathlib
import yaml
import time
import pyBigWig
from pybedtools import BedTool
from concurrent.futures import ProcessPoolExecutor, as_completed
from .utilities import determine_engine
from .region_ds_utilities import calculate_chunk_regions
import os
from ALLCools.utilities import parse_chrom_size
from .region_ds_utilities import update_region_ds_config

os.environ['NUMEXPR_MAX_THREADS'] = '16'


def _bigwig_over_bed(bed, path, value_type='mean', dtype='float32'):
    with pyBigWig.open(path, 'r') as bw:
        def _region_stat(row, t=value_type):
            chrom, start, end = row
            value = bw.stats(chrom, start, end, type=t)[0]
            return value

        values = bed.apply(_region_stat, t=value_type, axis=1)
        values = values.astype(dtype)
    return values


def _region_bed_sorted(bed_path, g, bed_sorted):
    chrom_sizes = parse_chrom_size(g)

    bed_df = pd.read_csv(bed_path, sep='\t', index_col=None, header=None)
    # select chroms that exist in g
    bed_df = bed_df.loc[bed_df.iloc[:, 0].isin(chrom_sizes.keys())]
    bed = BedTool.from_dataframe(bed_df)

    if bed_sorted:
        return bed
    else:
        return bed.sort(g=g)


def _bed_intersection(bed, path, g, region_index, bed_sorted):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        query_bed = _region_bed_sorted(path, g, bed_sorted)
        try:
            df = bed.intersect(query_bed,
                               wa=True,
                               f=0.2,
                               g=g,
                               sorted=True).to_dataframe()
            if df.shape[0] == 0:
                regions_idx = pd.Series([])
            else:
                regions_idx = df['name']
        except pd.errors.EmptyDataError:
            regions_idx = pd.Series([])
    regions = pd.Index(regions_idx.values)
    bool_series = pd.Series(region_index.isin(regions), index=region_index)

    query_bed.delete_temporary_history(ask=False)
    return bool_series


def _annotate_by_bigwigs_worker(
        dataset_path,
        region_dim,
        use_regions,
        chrom_size_path,
        bigwig_table,
        output_path,
        dim,
        slop,
        value_type,
        dtype, **kwargs):
    len(kwargs)
    # set dask scheduler to allow multiprocessing
    with dask.config.set(scheduler='sync'):
        # Open region ds again inside the worker function
        region_ds = RegionDS.open(path=dataset_path,
                                  region_dim=region_dim,
                                  use_regions=use_regions,
                                  chrom_size_path=chrom_size_path,
                                  select_dir=[region_dim])

        # get dmr region bed and bigwig files
        dmr_bed = region_ds.get_bed(with_id=False,
                                    bedtools=False,
                                    slop=slop,
                                    chrom_size_path=chrom_size_path)
        bigwig_files = pd.read_csv(bigwig_table,
                                   header=None,
                                   index_col=0,
                                   squeeze=True)

        # iterate each bigwig
        total_values = {}
        for sample, bigwig_path in bigwig_files.items():
            values = _bigwig_over_bed(bed=dmr_bed,
                                      path=bigwig_path,
                                      value_type=value_type,
                                      dtype=dtype)
            total_values[sample] = values
        total_values = pd.DataFrame(total_values)
        total_values.columns.name = dim

        ds = xr.Dataset({f'{region_dim}_{dim}_da': total_values}).sel({region_dim: use_regions})
        ds.to_zarr(output_path)
    return output_path


def _annotate_by_beds_worker(dataset_path, region_dim, use_regions, chrom_size_path,
                             slop, bed_table, dtype, dim, output_path, bed_sorted, **kwargs):
    len(kwargs)

    # set dask scheduler to allow multiprocessing
    with dask.config.set(scheduler='sync'):
        # Open region ds again inside the worker function
        region_ds = RegionDS.open(path=dataset_path,
                                  region_dim=region_dim,
                                  use_regions=use_regions,
                                  chrom_size_path=chrom_size_path,
                                  select_dir=[region_dim])

        # get dmr region bed
        dmr_bed = region_ds.get_bed(with_id=True,
                                    bedtools=True,
                                    slop=slop,
                                    chrom_size_path=chrom_size_path).sort(g=chrom_size_path)

        bed_files = pd.read_csv(bed_table,
                                header=None,
                                index_col=0,
                                squeeze=True)

        total_values = {}
        for sample, bed_path in bed_files.items():
            values = _bed_intersection(bed=dmr_bed,
                                       path=bed_path,
                                       bed_sorted=bed_sorted,
                                       g=chrom_size_path,
                                       region_index=region_ds.get_index(
                                           region_ds.get_region_dim()))
            total_values[sample] = values.astype(dtype)
        total_values = pd.DataFrame(total_values)
        total_values.columns.name = dim

        ds = xr.Dataset({f'{region_dim}_{dim}_da': total_values}).sel({region_dim: use_regions})
        ds.to_zarr(output_path)

        dmr_bed.delete_temporary_history(ask=False)
    return output_path


def _fisher_exact(row, alternative='two-sided'):
    from scipy.stats import fisher_exact
    oddsratio, p = fisher_exact(row.values.reshape((2, 2)), alternative=alternative)
    value = pd.Series({'oddsratio': oddsratio, 'p': p})
    return value


class RegionDS(xr.Dataset):
    __slots__ = ()

    def __init__(self, dataset, region_dim=None):
        super().__init__(data_vars=dataset.data_vars,
                         coords=dataset.coords,
                         attrs=dataset.attrs)
        if region_dim is not None:
            dataset.attrs['region_dim'] = region_dim
        return

    def set_region_dim(self, region_dim):
        if region_dim not in self.dims:
            raise KeyError(f'{region_dim} does not occur in dimension names: {list(self.dims.keys())}')
        self.attrs['region_dim'] = region_dim
        return

    def get_region_dim(self):
        return self.attrs['region_dim']

    def get_location(self):
        return self.attrs.get('region_ds_location')

    @classmethod
    def open(cls,
             path,
             region_dim=None,
             use_regions=None,
             split_large_chunks=True,
             chrom_size_path=None,
             select_dir=None):
        _path = pathlib.Path(path).absolute()
        region_ds_location = None
        region_ds_config = None

        if _path.is_dir():
            # check if this is a RegionDS dir that contains multiple datasets
            region_ds_file = _path / '.ALLCools'
            zarr_attr_file = _path / '.zattrs'
            exclude_dir_name = []
            if region_ds_file.exists():
                with open(region_ds_file) as f:
                    region_ds_config = yaml.load(f, yaml.SafeLoader)
                    if region_dim is None:
                        region_dim = region_ds_config['region_dim']
                        print(f'Using {region_dim} as region_dim')
                    # only open datasets having the region_dim
                    # other datasets will not be opened
                    # e.g, when region_dim == 'dmr', dms dataset will not be opened.
                    exclude_dir_name = [k for k, v in region_ds_config['ds_region_dim'].items() if v != region_dim]

                # add output_dir
                region_ds_location = pathlib.Path(path).absolute()

                # add chrom size path if exist
                if chrom_size_path is None:
                    chrom_size_path = _path / 'chrom_sizes.txt'
                    if not chrom_size_path.exists():
                        chrom_size_path = None

            if region_ds_file.exists() and not zarr_attr_file.exists():
                is_multi_region_ds = True
            else:
                is_multi_region_ds = False

            if is_multi_region_ds:
                # read all sub dir as a RegionDS, then merge all the RegionDS together.
                # dimension and coords should be compatible if all created by ALLCools
                region_dim_da = cls._open_single_dataset(path=f'{_path}/{region_dim}/*.zarr',
                                                         region_dim=region_dim,
                                                         split_large_chunks=split_large_chunks)
                datasets = [region_dim_da]
                for sub_dir_path in _path.iterdir():
                    if sub_dir_path.is_dir() and sub_dir_path.name[0] != '.':
                        if sub_dir_path.name == region_dim:
                            # already opened, skip
                            continue

                        if select_dir is not None and sub_dir_path.name not in select_dir:
                            # not in select_dir list, skip
                            continue

                        if sub_dir_path.name in exclude_dir_name:
                            # the dir do not have region_dim, skip
                            continue

                        if len(list(sub_dir_path.glob('*.zarr'))) == 0:
                            print(f'{sub_dir_path} does not seem to have any dataset, skipped')
                            continue
                        try:
                            datasets.append(
                                cls._open_single_dataset(path=f'{sub_dir_path}/*.zarr',
                                                         region_dim=region_dim,
                                                         split_large_chunks=split_large_chunks)
                            )
                        except BaseException as e:
                            print(f'An error raised when reading {sub_dir_path}.')
                            raise e
                total_dataset = xr.merge(datasets)
                total_dataset.attrs['region_dim'] = region_dim
                if chrom_size_path is not None:
                    total_dataset.attrs['chrom_size_path'] = chrom_size_path
                if region_ds_location is not None:
                    total_dataset.attrs['region_ds_location'] = region_ds_location
                region_ds = cls(total_dataset).squeeze()
            else:
                # is dir, but is not RegionDS dir, could be just a zarr dir
                region_ds = cls._open_single_dataset(path=path,
                                                     region_dim=region_dim,
                                                     split_large_chunks=split_large_chunks,
                                                     chrom_size_path=chrom_size_path)
        else:
            region_ds = cls._open_single_dataset(path=path,
                                                 region_dim=region_dim,
                                                 split_large_chunks=split_large_chunks,
                                                 chrom_size_path=chrom_size_path)

        if use_regions is not None:
            region_ds = region_ds.sel({region_dim: use_regions})

        # add config info into attrs
        if region_ds_config is not None:
            region_ds.attrs.update(region_ds_config)
        return region_ds

    @classmethod
    def _open_single_dataset(cls, path, region_dim, split_large_chunks=True, chrom_size_path=None):
        """
        Take one or multiple RegionDS file paths and create single RegionDS concatenated on region_dim

        Parameters
        ----------
        path
            Single RegionDS path or RegionDS path pattern with wildcard or RegionDS path list
        region_dim
            Dimension name of regions
        split_large_chunks
            Split large dask array chunks if true
        Returns
        -------
        RegionDS
        """
        # print('opening', path)
        if region_dim is None:
            raise ValueError('Please specify a region_dim name when open a normal xr.Dataset with RegionDS.')

        engine = determine_engine(path)
        # if engine is None:
        #     print(f'Open RegionDS with netcdf4 engine.')
        # else:
        #     print(f'Open RegionDS with {engine} engine.')
        if isinstance(path, str) and '*' not in path:
            ds = xr.open_dataset(path, engine=engine)
        else:
            with dask.config.set(
                    **{'array.slicing.split_large_chunks': split_large_chunks
                       }):
                if isinstance(path, str):
                    import glob
                    path = sorted([p for p in glob.glob(path)])
                ds = xr.open_mfdataset(path,
                                       parallel=False,
                                       combine='nested',
                                       concat_dim=region_dim,
                                       engine=engine)
        ds.attrs['region_dim'] = region_dim
        if chrom_size_path is not None:
            ds.attrs['chrom_size_path'] = chrom_size_path
        return cls(ds).squeeze()

    def iter_index(self, chunk_size=100000, dim=None):
        if dim is None:
            dim = self.get_region_dim()

        index = self.get_index(dim)
        for chunk_start in range(0, index.size, chunk_size):
            use_index = index[chunk_start: chunk_start + chunk_size]
            yield use_index

    def iter_array(self, chunk_size=100000, dim=None, da=None, load=False):
        if dim is None:
            dim = self.get_region_dim()

        if da is None:
            da = f'{dim}_da'
        _da = self[da]

        for _index in self.iter_index(chunk_size=chunk_size, dim=dim):
            use_da = _da.sel({dim: _index})
            if load:
                use_da.load()
            yield use_da

    def get_fasta(self,
                  genome_fasta,
                  output_path,
                  slop=None,
                  chrom_size_path=None,
                  standardize_length=None):
        bed = self.get_bed(with_id=True,
                           bedtools=True,
                           slop=slop,
                           chrom_size_path=chrom_size_path,
                           standardize_length=standardize_length)
        bed.getfasta(fo=output_path, nameOnly=True, fi=genome_fasta)
        return

    def get_bed(self,
                with_id=True,
                bedtools=False,
                slop=None,
                chrom_size_path=None,
                standardize_length=None):
        if chrom_size_path is None:
            chrom_size_path = self.attrs.get('chrom_size_path')  # will be none if not exist

        region_dim = self.get_region_dim()

        bed_df = pd.DataFrame({
            'chrom': self.coords[f'{region_dim}_chrom'],
            'start': self.coords[f'{region_dim}_start'],
            'end': self.coords[f'{region_dim}_end']
        })

        # standardize region length, used in motif enrichment analysis
        if standardize_length is not None:
            # standardize_length is an int number
            region_center = bed_df['start'] + (bed_df['end'] - bed_df['start']) // 2
            bed_df['start'] = region_center - 1
            bed_df['end'] = region_center
            slop = standardize_length // 2  # use the bedtools slop to extend the center to standard length

        if with_id:
            bed_df['name'] = self.get_index(region_dim).tolist()

        bed = None
        if slop is not None:
            if chrom_size_path is None:
                raise ValueError(
                    'Must provide chrom_size_path when slop is not None.')
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                bed = BedTool.from_dataframe(bed_df).slop(b=slop,
                                                          g=chrom_size_path)
            if not bedtools:
                bed_df = bed.to_dataframe()

        if bedtools:
            if bed is None:
                bed = BedTool.from_dataframe(bed_df)
            return bed
        else:
            bed_df.index = self.get_index(self.get_region_dim())
            return bed_df

    def _chunk_annotation_executor(self, annotation_function, cpu, **kwargs):
        chrom_size_path = kwargs['chrom_size_path']
        dim = kwargs['dim']
        try:
            input_table = kwargs['bigwig_table']
        except KeyError:
            input_table = kwargs['bed_table']
        chunk_size = kwargs['chunk_size']
        chunk_gbs = kwargs['chunk_gbs']
        dtype = kwargs['dtype']

        region_ds_path = self.get_location()
        if region_ds_path is None:
            raise ValueError(f'Must have an on-disk location to annotate bigwigs.')

        if chrom_size_path is None:
            chrom_size_path = self.attrs.get('chrom_size_path')
        region_dim = self.get_region_dim()
        n_regions = self.dims[region_dim]

        # deal with path
        dataset_path = f'{region_ds_path}/{region_dim}/*.zarr'
        output_dir = f'{region_ds_path}/{region_dim}_{dim}'
        _output_dir = pathlib.Path(output_dir)
        if _output_dir.exists():
            subprocess.run(f'rm -rf {output_dir}', shell=True)
        _output_dir.mkdir(exist_ok=True)

        n_features = pd.read_csv(input_table, header=None).shape[0]
        if chunk_size == 'auto':
            chunk_size = calculate_chunk_regions(chunk_size_gbs=chunk_gbs,
                                                 dtype=dtype,
                                                 n_features=n_features)
        print(f'Use chunk size {chunk_size}')
        n_chunks = n_regions // chunk_size + 1

        other_kwargs = {
            'region_dim': region_dim,
            'dataset_path': dataset_path,
            'chrom_size_path': chrom_size_path
        }
        kwargs.update(other_kwargs)

        with ProcessPoolExecutor(cpu) as exe:
            futures = {}
            region_chunk_iter = self.iter_index(chunk_size=chunk_size, dim=region_dim)
            for i, region_chunk in enumerate(region_chunk_iter):
                output_path = f'{output_dir}/chunk_{i:05d}.zarr'
                kwargs['output_path'] = output_path
                kwargs['use_regions'] = region_chunk
                future = exe.submit(annotation_function,
                                    **kwargs)
                futures[future] = i
                time.sleep(1)

            for i, future in enumerate(as_completed(futures)):
                chunk_i = futures[future]
                print(f'Job {i + 1}/{n_chunks} returned')
                output_path = future.result()
                print(f'Chunk {chunk_i} saved in {output_path}')

        # update config
        update_region_ds_config(output_dir, new_dataset_dim={dim: region_dim}, change_region_dim=None)

        # finally read the dataset into region ds
        _ds = RegionDS.open(f'{output_dir}/*.zarr', region_dim=region_dim)
        self.update(_ds)
        return

    def annotate_by_bigwigs(self,
                            bigwig_table,
                            dim,
                            slop=100,
                            chrom_size_path=None,
                            value_type='mean',
                            chunk_size='auto',
                            chunk_gbs=2,
                            dtype='float32',
                            cpu=1):
        kwargs = dict(
            bigwig_table=bigwig_table,
            dim=dim,
            slop=slop,
            chrom_size_path=chrom_size_path,
            value_type=value_type,
            chunk_size=chunk_size,
            chunk_gbs=chunk_gbs,
            dtype=dtype
        )
        self._chunk_annotation_executor(_annotate_by_bigwigs_worker, cpu, **kwargs)
        return

    def annotate_by_bed(self,
                        bed_table,
                        dim,
                        slop=100,
                        chrom_size_path=None,
                        chunk_size='auto',
                        chunk_gbs=4,
                        dtype='bool',
                        bed_sorted=True,
                        cpu=1):
        bed_tmp = pathlib.Path(f'./pybedtools_tmp_{np.random.randint(0, 100000)}').absolute()
        bed_tmp.mkdir(exist_ok=True)
        pybedtools.helpers.set_tempdir(bed_tmp)

        kwargs = dict(
            bed_table=bed_table,
            dim=dim,
            slop=slop,
            chrom_size_path=chrom_size_path,
            chunk_size=chunk_size,
            chunk_gbs=chunk_gbs,
            dtype=dtype,
            bed_sorted=bed_sorted
        )
        self._chunk_annotation_executor(_annotate_by_beds_worker, cpu, **kwargs)

        subprocess.run(f'rm -rf {bed_tmp}', shell=True)
        return

    def scan_motifs(self,
                    genome_fasta,
                    cpu=1,
                    standardize_length=500,
                    motif_set_path=None,
                    chrom_size_path=None,
                    combine_cluster=True,
                    fnr_fpr_fold=1000,
                    chunk_size=None,
                    dim='motif'):
        region_ds_path = self.get_location()
        if region_ds_path is None:
            raise ValueError(f'Must have an on-disk location to annotate bigwigs.')

        if chrom_size_path is None:
            chrom_size_path = self.attrs.get('chrom_size_path')
        region_dim = self.get_region_dim()

        fasta_path = f'{region_ds_path}/regions.fasta'
        self.get_fasta(genome_fasta,
                       output_path=fasta_path,
                       slop=None,
                       chrom_size_path=chrom_size_path,
                       standardize_length=standardize_length)

        from ..motif import MotifSet
        if motif_set_path is not None:
            motif_set: MotifSet = joblib.load(motif_set_path)
        else:
            # default motif database from three sources
            from ..motif import get_default_motif_set
            motif_set = get_default_motif_set()

        if fnr_fpr_fold != 1000:
            motif_set.calculate_threshold(cpu=cpu, method='balance', threshold_value=fnr_fpr_fold)

        motif_ds = motif_set.scan_motifs(
            fasta_path=fasta_path,
            output_dir=region_ds_path,
            cpu=cpu,
            region_dim=region_dim,
            combine_cluster=combine_cluster,
            motif_dim=dim,
            chunk_size=chunk_size)

        self.update(motif_ds)
        return

    def get_hypo_hyper_index(self,
                             a,
                             region_dim=None,
                             region_state_da=None,
                             sample_dim='sample'):
        if region_state_da is None:
            if region_dim is None:
                region_dim = self.get_region_dim()
            region_state_da = f'{region_dim}_state'
        a_states = self[region_state_da].sel({sample_dim: a}).to_pandas()

        hypo_dmr = a_states == -1
        hypo_dmr = hypo_dmr[hypo_dmr].index

        hyper_dmr = a_states == 1
        hyper_dmr = hyper_dmr[hyper_dmr].index

        return hypo_dmr, hyper_dmr

    def get_pairwise_differential_index(self,
                                        a,
                                        b,
                                        dmr_type='hypo',
                                        region_dim=None,
                                        region_state_da=None,
                                        sample_dim='sample'):
        a_hypo, a_hyper = self.get_hypo_hyper_index(a,
                                                    region_dim=region_dim,
                                                    region_state_da=region_state_da,
                                                    sample_dim=sample_dim)
        b_hypo, b_hyper = self.get_hypo_hyper_index(b,
                                                    region_dim=region_dim,
                                                    region_state_da=region_state_da,
                                                    sample_dim=sample_dim)

        if dmr_type.lower() == 'hypo' or dmr_type == -1:
            # a hypo, b not hypo
            a_not_b = a_hypo[~a_hypo.isin(b_hypo)]
            a_and_b = a_hypo[a_hypo.isin(b_hypo)]
            b_not_a = b_hypo[~b_hypo.isin(a_hypo)]
        elif dmr_type.lower() == 'hyper' or dmr_type == 1:
            # a hypo, b not hypo
            a_not_b = a_hyper[~a_hyper.isin(b_hyper)]
            a_and_b = a_hyper[a_hyper.isin(a_hyper)]
            b_not_a = b_hyper[~b_hyper.isin(a_hyper)]
        else:
            raise ValueError(f'Unknown dmr_type {dmr_type}.')

        return a_not_b, a_and_b, b_not_a

    def motif_enrichment(self,
                         true_regions,
                         background_regions,
                         region_dim=None,
                         motif_dim='motif-cluster',
                         motif_da=None,
                         alternative='two-sided'):
        if region_dim is None:
            region_dim = self.get_region_dim()
        if motif_da is None:
            motif_da = f'{region_dim}_{motif_dim}_da'
        true_motif = self[motif_da].sel({'motif_value': 'n_motifs',
                                         region_dim: true_regions})
        true_motif = (true_motif > 0).sum(dim=region_dim).to_pandas()
        true_no_motif = true_regions.size - true_motif

        bkg_motif = self[motif_da].sel({'motif_value': 'n_motifs',
                                        region_dim: background_regions})
        bkg_motif = (bkg_motif > 0).sum(dim=region_dim).to_pandas()
        bkg_no_motif = background_regions.size - bkg_motif

        contingency_tables = pd.DataFrame({
            'tp': true_motif,
            'tn': true_no_motif,
            'fp': bkg_motif,
            'fn': bkg_no_motif
        })

        # add pseudo count
        contingency_tables += 1

        test_results = contingency_tables.apply(_fisher_exact, axis=1, alternative=alternative)
        from statsmodels.stats.multitest import multipletests
        _, q, *_ = multipletests(test_results['p'], method='fdr_bh')
        test_results['q'] = q

        test_results['log2OR'] = np.log2(test_results['oddsratio'])
        test_results['-lgq'] = -np.log10(q)

        return test_results

    def sample_dmr_motif_enrichment(self,
                                    sample,
                                    region_dim=None,
                                    sample_dim='sample',
                                    motif_dim='motif-cluster',
                                    region_state_da=None,
                                    motif_da=None,
                                    alternative='two-sided'):
        hypo_region, hyper_region = self.get_hypo_hyper_index(sample,
                                                              region_dim=region_dim,
                                                              region_state_da=region_state_da,
                                                              sample_dim=sample_dim)
        test_results = self.motif_enrichment(hypo_region,
                                             hyper_region,
                                             region_dim=region_dim,
                                             motif_dim=motif_dim,
                                             motif_da=motif_da,
                                             alternative=alternative)
        return test_results

    def pairwise_dmr_motif_enrichment(self, a, b, dmr_type='hypo', region_dim=None, region_state_da=None,
                                      sample_dim='sample', motif_dim='motif-cluster', motif_da=None,
                                      alternative='two-sided'):
        a_not_b, _, b_not_a = self.get_pairwise_differential_index(a,
                                                                   b,
                                                                   dmr_type=dmr_type,
                                                                   region_dim=region_dim,
                                                                   region_state_da=region_state_da,
                                                                   sample_dim=sample_dim)
        test_results = self.motif_enrichment(a_not_b,
                                             b_not_a,
                                             region_dim=region_dim,
                                             motif_dim=motif_dim,
                                             motif_da=motif_da,
                                             alternative=alternative)
        return test_results

    def region_correlation(self):
        return

    def object_coords_to_string(self):
        for k, v in self.coords.items():
            if v.dtype == 'object':
                self.coords[k] = v.load().astype(str)
        return

    def save(self, path):
        # turn object coords to fix length string dtype before saving to zarr
        self.object_coords_to_string()

        if pathlib.Path(path).exists():
            subprocess.run(f'rm -rf {path}', shell=True)
        self.to_zarr(path)

        with open(f'{path}/.ALLCools', 'w') as f:
            yaml.dump(self.attrs, f)
        return

    def get_coords(self, name):
        return self.coords[name].to_pandas()

import pysam
import pandas as pd
import numpy as np
from scipy import stats
from scipy.sparse import csr_matrix
import pathlib
import anndata
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from ..utilities import parse_mc_pattern


def _read_region_bed(bed_path):
    region_bed = pd.read_csv(bed_path, sep='\t', header=None, index_col=3)
    region_bed.index.name = 'region'
    region_bed.columns = ['chrom', 'start', 'end']
    return region_bed


def _count_single_allc(allc_path, bed_path, mc_pattern, output_dir):
    patterns = parse_mc_pattern(mc_pattern)
    region_bed = _read_region_bed(bed_path)

    # bin raw counts
    with pysam.TabixFile(allc_path, 'r') as allc:
        records = []  # list of [mc, cov, region_idx]
        for idx, (region_id, (chrom, start,
                              end)) in enumerate(region_bed.iterrows()):
            total_mc = 0
            total_cov = 0
            try:
                iterator = allc.fetch(chrom, start, end)
            except ValueError:
                # in low coverage cells, allc file might miss a whole chromosome,
                # this cause value error saying "could not create iterator for region"
                continue
            for line in iterator:
                chrom, pos, _, context, mc, cov, _ = line.split('\t')
                if context in patterns:
                    total_mc += int(mc)
                    total_cov += int(cov)
            if total_cov > 0:
                records.append([idx, total_mc, total_cov])
        bin_counts = pd.DataFrame(records, columns=['idx', 'mc',
                                                    'cov']).set_index('idx')

    # calculate binom sf (1-cdf) value, hypo bins are close to 1, hyper bins are close to 0
    mc_sum, cov_sum = bin_counts.sum()
    p = mc_sum / (cov_sum + 0.000001)  # prevent empty allc error
    pv = bin_counts.apply(lambda x: stats.binom(x['cov'], p).sf(x['mc']) if x['mc'] < x['cov'] else 0,
                          axis=1).astype('float16')
    pv = pv[pv > 0.5]  # get rid of most hyper bins
    pv.to_hdf(f'{output_dir}/{pathlib.Path(allc_path).name}.hdf', key='data')
    return


def generate_mcad(allc_table, bed_path, output_path, mc_pattern, cpu, cleanup=True):
    """
    Generate cell-by-region sparse matrix for a single methylation pattern.

    Parameters
    ----------
    allc_table
    bed_path
    output_path
    mc_pattern
    cpu
    cleanup

    Returns
    -------

    """
    # allc table has 2 columns: cell_id \t allc_path
    allc_paths = pd.read_csv(allc_table,
                             sep='\t',
                             index_col=0,
                             header=None,
                             squeeze=True)
    allc_paths.index.name = 'cell'

    # temp dir
    temp_dir = pathlib.Path(f'{output_path}_pv_temp')
    temp_dir.mkdir(exist_ok=True)

    # calculating individual cells
    with ProcessPoolExecutor(cpu) as executor:
        futures = {}
        allc_path_idy = {}
        for idy, (cell_id,
                  allc_path) in enumerate(allc_paths.items()):
            allc_path_idy[pathlib.Path(allc_path).name] = idy
            output_path = temp_dir / f'{pathlib.Path(allc_path).name}.hdf'
            if output_path.exists():
                continue
            future = executor.submit(_count_single_allc,
                                     allc_path=allc_path,
                                     bed_path=bed_path,
                                     mc_pattern=mc_pattern,
                                     output_dir=temp_dir)
            futures[future] = cell_id

        for future in as_completed(futures):
            cell_id = futures[future]
            print(f'{cell_id} returned.')
            future.result()

    # aggregate all the cells
    print('Aggregate cells into adata')
    total_idx = []
    total_idy = []
    total_data = []
    for path in temp_dir.glob('*hdf'):
        idy = allc_path_idy[path.name[:-4]]
        pv = pd.read_hdf(path)
        if pv.size == 0:
            # no sig result
            continue
        total_idx.append(pv.index.values)
        total_idy.append([idy] * pv.size)
        total_data.append(pv.values)
    # the cell by region matrix
    region_bed = _read_region_bed(bed_path)
    _data = csr_matrix(
        (np.concatenate(total_data),
         (np.concatenate(total_idy), np.concatenate(total_idx))),
        shape=(len(allc_path_idy), region_bed.shape[0]))

    # save the data as anndata
    adata = anndata.AnnData(_data,
                            obs=pd.DataFrame([], index=allc_paths.index),
                            var=region_bed)
    adata.write_h5ad(f'{output_path}.mcad')

    # remove temp
    if cleanup:
        subprocess.run(['rm', '-rf', str(temp_dir)])
    return

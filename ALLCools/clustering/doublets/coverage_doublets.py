import pandas as pd
from collections import defaultdict
import pysam
import pathlib
import json
import numpy as np
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
from concurrent.futures import ProcessPoolExecutor, as_completed


def _calculate_cell_record(allc_path, output_path, cov_cutoff=2, resolution=100):
    """Count the high coverage bins for each cell, save results to json"""
    allc_path = str(allc_path)
    output_path = str(output_path)
    cov_high_cutoff = int(cov_cutoff * 2)

    cell_records = {}
    i = 0
    with pysam.TabixFile(allc_path) as allc:
        for i, line in enumerate(allc.fetch()):
            chrom, pos, *_, cov, _ = line.split('\t')
            cov = int(cov)
            if cov_cutoff < cov <= cov_high_cutoff:
                bin_id = int(pos) // resolution
                try:
                    cell_records[chrom][bin_id] += 1
                except KeyError:
                    cell_records[chrom] = defaultdict(int)
                    cell_records[chrom][bin_id] += 1

    cell_total_c = i + 1
    cell_records = {
        chrom: list(values.keys())
        for chrom, values in cell_records.items()
    }
    # final output to disk
    total_records = {'total_c': cell_total_c, 'bins': cell_records}

    with open(output_path, 'w') as f:
        json.dump(total_records, f)
    return output_path


def calculate_blacklist_region(region_records, alpha=0.01):
    """Collect highly covered regions by region-wise poisson FDR p value < alpha"""
    # calculate region poisson mu
    sum_of_bin = 0
    n_bin = 0
    for chrom, chrom_values in region_records.items():
        sum_of_bin += sum(chrom_values.values())
        n_bin += len(chrom_values)
    mu = sum_of_bin / n_bin

    # calculate region FDR p cutoff
    total_p = []
    for chrom, chrom_values in region_records.items():
        chrom_values = pd.Series(chrom_values)
        p_values = poisson.sf(chrom_values.values, mu)
        total_p.append(p_values)
    total_p = np.concatenate(total_p)
    judge, *_ = multipletests(total_p, alpha=alpha, method='fdr_bh')
    p_max = total_p[judge].max()
    del total_p, judge

    # calculate region blacklist
    final_blacklist = {}
    for chrom, chrom_values in region_records.items():
        chrom_values = pd.Series(chrom_values)
        p_values = poisson.sf(chrom_values.values, mu)
        final_blacklist[chrom] = list(chrom_values[p_values < p_max].index)
    return final_blacklist


def _calculate_cell_final_values(output_path, region_blacklist):
    """Calculate final cell values while remove blacklist"""
    with open(output_path) as f:
        cell_record = json.load(f)
        total_n = 0
        for chrom, bins in cell_record['bins'].items():
            total_n += len(set(bins) - region_blacklist[chrom])
    return total_n, cell_record['total_c']


def coverage_doublets(allc_dict: dict,
                      resolution: int = 100,
                      cov_cutoff=2,
                      region_alpha=0.01,
                      tmp_dir='doublets_temp_dir',
                      cpu=1,
                      keep_tmp=False):
    """
    Quantify cell high coverage bins for doublets evaluation

    Parameters
    ----------
    allc_dict
        dict with cell_id as key, allc_path as value
    resolution
        genome bin resolution to quantify, bps
    cov_cutoff
        cutoff the cov, sites within cov_cutoff < cov <= 2 * cov_cutoff will be count
    region_alpha
        FDR adjusted P-value cutoff
    tmp_dir
        temporary dir to save the results
    cpu
        number of cpu to use
    keep_tmp
        Whether save the tem_dir for debugging

    Returns
    -------

    """
    tmp_dir = pathlib.Path(tmp_dir)
    tmp_dir.mkdir(exist_ok=True)

    # count each cell and collect region-wise sum in the same time
    region_records = {}

    def _sum_region(p):
        with open(p) as cr:
            cell_record = json.load(cr)
            for chrom, chrom_bins in cell_record['bins'].items():
                if chrom not in region_records:
                    region_records[chrom] = defaultdict(int)
                for bin_id in chrom_bins:
                    region_records[chrom][bin_id] += 1
        return

    cell_paths = {}
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for cell_id, path in allc_dict.items():
            output_path = f"{tmp_dir}/{cell_id}.json"
            if pathlib.Path(output_path).exists():
                # directly quantify region records
                _sum_region(output_path)
                cell_paths[cell_id] = output_path
                continue
            future = exe.submit(_calculate_cell_record,
                                allc_path=path,
                                output_path=output_path,
                                resolution=resolution,
                                cov_cutoff=cov_cutoff)
            futures[future] = cell_id

        # during calculating the cell records, also summarize region records
        for future in as_completed(futures):
            cell_id = futures[future]
            output_path = future.result()
            cell_paths[cell_id] = output_path
            _sum_region(output_path)

    # calculate dataset specific region blacklist
    region_blacklist = calculate_blacklist_region(region_records, alpha=region_alpha)
    with open(f'{tmp_dir}/region_blacklist.json', 'w') as f:
        json.dump(region_blacklist, f)
    # list to set, dump don't support set
    region_blacklist = {k: set(v) for k, v in region_blacklist.items()}

    # calculate cell final stats
    total_values = []
    cells = []
    for cell_id, output_path in cell_paths.items():
        cell_values = _calculate_cell_final_values(output_path, region_blacklist)
        total_values.append(cell_values)
        cells.append(cell_id)
    total_values = pd.DataFrame(total_values, index=cells, columns=['TotalHCB', 'TotalC'])
    # this value don't have specific mathematical meaning, for plotting
    total_values['HCBRatio'] = total_values['TotalHCB'] * resolution / total_values['TotalC']

    if not keep_tmp:
        import os
        os.rmdir(tmp_dir)

    return total_values

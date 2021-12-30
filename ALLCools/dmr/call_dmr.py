import xarray as xr
import numpy as np
import pandas as pd
from ..mcds import RegionDS
from ..mcds.region_ds_utilities import update_region_ds_config
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def call_dmr_single_chrom(output_dir,
                          chrom,
                          p_value_cutoff=0.001,
                          frac_delta_cutoff=0.2,
                          max_dist=250,
                          residual_quantile=0.6,
                          corr_cutoff=0.3):
    ds = RegionDS.open(f'{output_dir}/dms/*.zarr', region_dim='dms')

    # open DMS dataset and select single chromosome, significant, large delta DMS
    p_value_judge = ds.coords['dms_p-values'] < p_value_cutoff
    frac_delta = ds['dms_da_frac'].max(dim='sample') - ds['dms_da_frac'].min(dim='sample')
    frac_delta_judge = frac_delta > frac_delta_cutoff
    chrom_judge = ds['dms_chrom'] == chrom
    # has to use load(scheduler='sync') when this function is called from multiprocess
    final_judge = (p_value_judge & frac_delta_judge & chrom_judge).load(scheduler='sync').to_pandas()
    if final_judge.sum() == 0:
        # No DMS remain
        return

    use_dms = final_judge[final_judge].index
    ds = ds.sel(dms=use_dms).load(scheduler='sync')

    # step 1: combine raw DMR windows based on distance
    dms_dist = ds['dms_pos'][1:].values - ds['dms_pos'][:-1].values
    dist_judge = dms_dist < max_dist
    cur_dmr = 0
    dmr_ids = [0]
    for i in dist_judge:
        if not i:
            cur_dmr += 1
        dmr_ids.append(cur_dmr)
    dmr_ids = pd.Series(dmr_ids, index=ds.get_index('dms'))

    # step 2: determine correlation between adjacent DMS
    data = ds['dms_da_frac'].transpose('dms', 'sample').to_pandas()
    a = data.iloc[:-1, :].reset_index(drop=True)
    b = data.iloc[1:, :].reset_index(drop=True)
    # index of corr means the correlation of that DMS with the previous DMS
    # regardless of the genome distance
    # fill na value with 1, tend to not to split DMR due to nan
    corrs = a.fillna(1).corrwith(b.fillna(1), axis=1, method='pearson')
    corrs.index = data.iloc[1:, :].index.copy()

    # step 3: recalculate DMR windows based on both distance and correlation
    raw_dmr_table = pd.DataFrame({'dmr': dmr_ids, 'corr': corrs})
    dmr_dict = {}
    cur_dmr_id = 0
    for dmr_id, sub_df in raw_dmr_table.groupby('dmr'):
        dmr_dict[sub_df.index[0]] = cur_dmr_id
        if sub_df.shape[0] > 1:
            for dms_id, corr in sub_df['corr'][1:].items():
                if corr > corr_cutoff:
                    dmr_dict[dms_id] = cur_dmr_id
                else:
                    cur_dmr_id += 1
                    dmr_dict[dms_id] = cur_dmr_id
        cur_dmr_id += 1
    dmr_ids = pd.Series(dmr_dict)
    dmr_ids.index.name = 'dms'
    ds.coords['dmr'] = dmr_ids

    # step 4: determine sample hypo or hyper in each DMS and DMR based on residual
    if residual_quantile < 0.5:
        residual_quantile = 1 - residual_quantile

    residual_low_cutoff, residual_high_cutoff = np.nanquantile(
        ds['dms_residual'], [1 - residual_quantile, residual_quantile]
    )

    lower_residual = ds['dms_residual'] < residual_low_cutoff
    higher_residual = ds['dms_residual'] > residual_high_cutoff
    # for each sample in each DMS
    # -1 means hypo methylation, 1 means hyper methylation, 0 mean no sig change
    dms_states = lower_residual.astype(int) * -1 + higher_residual.astype(int)
    # dmr state judge
    dmr_states = dms_states.groupby('dmr').sum()
    dmr_states = xr.where(dmr_states > 1, 1, dmr_states)
    dmr_states = xr.where(dmr_states < -1, -1, dmr_states)

    # step 5: prepare dmr counts and fractions
    dmr_da = ds['dms_da'].groupby('dmr').sum()
    dmr_frac = dmr_da.sel(count_type='mc') / dmr_da.sel(count_type='cov')
    dmr_ds = xr.Dataset({
        'dmr_da': dmr_da.astype(np.uint32),
        'dmr_state': dmr_states.astype(np.float32),
        'dmr_da_frac': dmr_frac.astype(np.int8)
    })

    # add n dms counts
    n_dms = dmr_ids.value_counts().sort_index()
    n_dms.index.name = 'dmr'
    dmr_ds.coords['ndms'] = n_dms
    # add genome coords
    dmr_ds.coords['chrom'] = pd.Series([chrom for _ in range(dmr_ds.get_index('dmr').size)],
                                       index=dmr_ds.get_index('dmr'))
    dmr_ds.coords['start'] = ds['dms_pos'].groupby(ds['dmr']).min().to_pandas() - 1
    dmr_ds.coords['end'] = ds['dms_pos'].groupby(ds['dmr']).max().to_pandas() + 1
    dmr_ds.coords['length'] = dmr_ds.coords['end'] - dmr_ds.coords['start']
    dmr_ds.coords['dmr'] = dmr_ds['chrom'].to_pandas().astype(str) + '-' + dmr_ds['dmr'].to_pandas().astype(str)

    # rename none dimensional coords to prevent collision when merge with other ds
    dmr_ds = dmr_ds.rename({k: f'dmr_{k}' for k in dmr_ds.coords.keys() if k not in dmr_ds.dims})

    dmr_ds.to_zarr(f'{output_dir}/dmr/{chrom}.zarr')
    return


def call_dmr(output_dir,
             p_value_cutoff=0.001,
             frac_delta_cutoff=0.2,
             max_dist=250,
             residual_quantile=0.6,
             corr_cutoff=0.3,
             cpu=1,
             chrom=None):
    subprocess.run(f'rm -rf {output_dir}/dmr*', shell=True)
    update_region_ds_config(output_dir=output_dir,
                            new_dataset_dim={'dmr': 'dmr'},
                            change_region_dim='dmr')

    if chrom is None:
        chrom_size_path = f'{output_dir}/chrom_sizes.txt'
        from ..utilities import parse_chrom_size
        chroms = parse_chrom_size(chrom_size_path).keys()
    else:
        if isinstance(chrom, list):
            chroms = chrom
        else:
            chroms = [chrom]

    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for chrom in chroms:
            future = exe.submit(call_dmr_single_chrom,
                                output_dir=output_dir,
                                chrom=chrom,
                                p_value_cutoff=p_value_cutoff,
                                frac_delta_cutoff=frac_delta_cutoff,
                                max_dist=max_dist,
                                residual_quantile=residual_quantile,
                                corr_cutoff=corr_cutoff)
            futures[future] = chrom

        for future in as_completed(futures):
            chrom = futures[future]
            print(f'{chrom} returned.')
            future.result()
    return

import numpy as np


def calculate_posterior_mc_rate(mc_da, cov_da, var_dim,
                                normalize_per_cell=True, clip_norm_value=10):
    # TODO calculate cell_a, cell_b separately
    # TODO add a parameter weighting var to adjust prior
    # so we can do post_rate only in a very small set of gene to prevent memory issue

    raw_rate = mc_da / cov_da
    cell_rate_mean = raw_rate.mean(dim=var_dim)  # this skip na
    cell_rate_var = raw_rate.var(dim=var_dim)  # this skip na

    # based on beta distribution mean, var
    # a / (a + b) = cell_rate_mean
    # a * b / ((a + b) ^ 2 * (a + b + 1)) = cell_rate_var
    # calculate alpha beta value for each cell
    cell_a = (1 - cell_rate_mean) * (cell_rate_mean ** 2) / cell_rate_var - cell_rate_mean
    cell_b = cell_a * (1 / cell_rate_mean - 1)

    # cell specific posterior rate
    post_rate = (mc_da + cell_a) / (cov_da + cell_a + cell_b)

    if normalize_per_cell:
        # there are two ways of normalizing per cell, by posterior or prior mean:
        # prior_mean = cell_a / (cell_a + cell_b)
        # posterior_mean = post_rate.mean(dim=var_dim)

        # Here I choose to use prior_mean to normalize cell,
        # therefore all cov == 0 features will have normalized rate == 1 in all cells.
        # i.e. 0 cov feature will provide no info
        prior_mean = cell_a / (cell_a + cell_b)
        post_rate = post_rate / prior_mean
        if clip_norm_value is not None:
            post_rate.values[np.where(post_rate.values > clip_norm_value)] = clip_norm_value
    return post_rate


def calculate_gch_rate(mcds, var_dim='chrom100k'):
    rate_da = mcds.sel(mc_type=['GCHN', 'HCHN']).add_mc_rate(dim=var_dim, da=f'{var_dim}_da',
                                                             normalize_per_cell=False, inplace=False)
    # (PCG - PCH) / (1 - PCH)
    real_gc_rate = (rate_da.sel(mc_type='GCHN') - rate_da.sel(mc_type='HCHN')) / (
            1 - rate_da.sel(mc_type='HCHN'))
    real_gc_rate = real_gc_rate.transpose('cell', var_dim).values
    real_gc_rate[real_gc_rate < 0] = 0

    # norm per cell
    cell_overall_count = mcds[f'{var_dim}_da'].sel(mc_type=['GCHN', 'HCHN']).sum(dim=var_dim)
    cell_overall_rate = cell_overall_count.sel(count_type='mc') / cell_overall_count.sel(count_type='cov')
    gchn = cell_overall_rate.sel(mc_type='GCHN')
    hchn = cell_overall_rate.sel(mc_type='HCHN')
    overall_gchn = (gchn - hchn) / (1 - hchn)
    real_gc_rate = real_gc_rate / overall_gchn.values[:, None]
    return real_gc_rate

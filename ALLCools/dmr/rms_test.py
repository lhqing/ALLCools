import numpy as np
from numba import njit
from scipy.stats import norm


@njit
def _get_e(table, n, m):
    row_sum = table.sum(axis=0).astype(np.float64)
    col_sum = table.sum(axis=1).astype(np.float64)
    e = np.dot(col_sum.reshape(n, 1), row_sum.reshape(1, 2)) / m
    return e


@njit
def _calculate_goodness_of_fit(table, n, m):
    e = _get_e(table, n, m)
    s = np.sqrt(((table - e) ** 2).sum() / 2 / n)
    return s


@njit
def _make_permutation_table(p, n, m):
    permuted_flat = np.zeros_like(p, dtype=np.int64)
    order = np.arange(0, permuted_flat.size)
    for _ in range(m):
        permuted_flat[np.random.choice(order)] += 1
    return permuted_flat.reshape(n, 2)


@njit
def calculate_residual(table):
    """Calculate residual of the table."""
    n = table.shape[0]
    m = table.sum()
    e = _get_e(table, n, m)
    residual = (table - e) / np.sqrt(
        np.multiply(
            e,
            (1 - e.sum(axis=1) / m).reshape(n, 1) * (1 - e.sum(axis=0) / m).reshape(1, 2),
        )
    )
    return residual


@njit
def permute_root_mean_square_test(table, n_permute=3000, min_pvalue=0.034):
    """
    Permute root-mean-square test.

    Parameters
    ----------
    table :
        N-by-2 table of methylated and unmethylated counts.
    n_permute :
        Number of permutations.
    min_pvalue :
        Minimum p-value to perform the test.
    """
    # calculate real goodness-of-fit s
    n = table.shape[0]
    m = table.sum()
    e = _get_e(table, n, m)
    real_s = _calculate_goodness_of_fit(table, n, m)

    # permutation
    p = e.flatten() / m
    greater_than_real = 1
    max_greater_value = n_permute * min_pvalue
    for i in range(n_permute):
        p_table = _make_permutation_table(p, n, m)
        # calculate permuted goodness of fit s'
        s = _calculate_goodness_of_fit(p_table, n, m)
        greater_than_real += int(s >= real_s)
        # break in advance if p-value can be significant
        if greater_than_real > max_greater_value:
            # return current p value
            return greater_than_real / (i + 2)

    p_value = greater_than_real / n_permute
    return p_value


@njit
def _get_read_s_and_permute_s(table, n_permute=10000, min_pvalue=0.034):
    # calculate real goodness-of-fit s
    n = table.shape[0]
    m = table.sum()
    e = _get_e(table, n, m)
    real_s = _calculate_goodness_of_fit(table, n, m)

    # permutation
    p = e.flatten() / m
    greater_than_real = 1
    max_greater_value = n_permute * min_pvalue
    ss = np.zeros(n_permute + 1, dtype=np.float64)
    ss[0] = real_s
    for i in range(n_permute):
        p_table = _make_permutation_table(p, n, m)
        # calculate permuted goodness of fit s'
        s = _calculate_goodness_of_fit(p_table, n, m)
        ss[i + 1] = s

        greater_than_real += int(s >= real_s)
        # break in advance if p-value can be significant
        if greater_than_real > max_greater_value:
            # return current p value
            return ss[: i + 1].copy()
    return ss


def permute_root_mean_square_test_and_estimate_p(table, n_permute=10000, min_pvalue=0.034):
    """
    Permute root-mean-square test

    Report estimated p-value by approximate null distribution as normal distribution.
    """
    ss = _get_read_s_and_permute_s(table=table, n_permute=n_permute, min_pvalue=min_pvalue)
    real_s = ss[0]
    ss = ss[1:]

    e, std = norm.fit(ss)
    if std == 0:
        p_value = 1
    else:
        p_value = 1 - norm(e, std).cdf(real_s)

    if p_value < 1e-14:  # np.finfo(np.float64).resolution: 1e-15
        p_value = 1e-14
    return p_value


@njit
def _downsample_sample_count(a, max_count):
    a = a.astype(np.float64)
    total = a.sum()
    p = a / total
    if total > max_count:
        b = max_count * p
    else:
        b = a
    b = b.astype(np.int64)
    return b


def downsample_table(table, max_row_count=50, max_total_count=3000):
    """Downsample high count rows to max_row_count."""
    n_sample = table.shape[0]
    max_row_count = min(max_total_count / n_sample, max_row_count)
    max_row_count = max(10, max_row_count)

    return np.apply_along_axis(_downsample_sample_count, axis=1, arr=table, max_count=max_row_count)


def init_rms_functions():
    """Initialize the numba JIT compilation before multiprocessing RMS functions."""
    table = np.array([[0, 1], [0, 1]])
    permute_root_mean_square_test(table)
    calculate_residual(table)
    downsample_table(table, max_row_count=10)
    return

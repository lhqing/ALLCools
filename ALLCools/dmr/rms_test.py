import numpy as np
from numba import njit


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
    n = table.shape[0]
    m = table.sum()
    e = _get_e(table, n, m)
    residual = (table - e) / np.sqrt(
        np.multiply(
            e,
            (1 - e.sum(axis=1) / m).reshape(n, 1)
            * (1 - e.sum(axis=0) / m).reshape(1, 2),
        )
    )
    return residual


@njit
def permute_root_mean_square_test(table, n_permute=3000, min_pvalue=0.034):
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


def downsample_table(table, max_row_count):
    """Downsample high count rows to max_row_count"""
    return np.apply_along_axis(
        _downsample_sample_count, axis=1, arr=table, max_count=max_row_count
    )

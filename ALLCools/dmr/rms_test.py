import numpy as np
from numba import njit


@njit(cache=True)
def calculate_goodness_of_fit(table):
    n = table.shape[0]
    m = table.sum()
    row_sum = table.sum(axis=0).astype(np.float64)
    col_sum = table.sum(axis=1).astype(np.float64)
    e = np.dot(col_sum.reshape(n, 1), row_sum.reshape(1, 2)) / m
    s = np.sqrt(((table - e) ** 2).sum() / 2 / n)
    return s


@njit(cache=True)
def make_permutation_table(p, m, n):
    permuted_flat = np.zeros_like(p, dtype=np.int64)
    order = np.arange(0, permuted_flat.size)
    for _ in range(m):
        permuted_flat[np.random.choice(order)] += 1
    return permuted_flat.reshape(n, 2)


@njit(cache=True)
def permute_root_mean_square_test(table, n_permute=3000, min_pvalue=0.034):
    # calculate real goodness-of-fit s
    n = table.shape[0]
    m = table.sum()
    row_sum = table.sum(axis=0).astype(np.float64)
    col_sum = table.sum(axis=1).astype(np.float64)
    e = np.dot(col_sum.reshape(n, 1), row_sum.reshape(1, 2)) / m
    real_s = calculate_goodness_of_fit(table)

    # permutation
    p = e.flatten() / m
    greater_than_real = 1
    max_greater_value = n_permute * min_pvalue
    for i in range(n_permute):
        p_table = make_permutation_table(p, m, n)
        # calculate permuted goodness of fit s'
        s = calculate_goodness_of_fit(p_table)
        greater_than_real += int(s >= real_s)
        # break in advance if p-value can be significant
        if greater_than_real > max_greater_value:
            # return current p value
            return greater_than_real / (i + 2)

    p_value = greater_than_real / n_permute
    return p_value


@njit
def calculate_residue(table):
    n = table.shape[0]
    m = table.sum()
    row_sum = table.sum(axis=0).astype(np.float64)
    col_sum = table.sum(axis=1).astype(np.float64)
    e = np.dot(col_sum.reshape(n, 1), row_sum.reshape(1, 2)) / m
    residual = (table - e) / np.sqrt(
        np.multiply(e, (1 - e.sum(axis=1) / m).reshape(n, 1) *
                    (e.sum(axis=0) / m).reshape(1, 2)))
    return residual

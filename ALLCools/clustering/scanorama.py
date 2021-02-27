"""
Allow to chose different metric

Modified from Scanorama, see LICENSE
https://github.com/brianhie/scanorama/blob/master/LICENSE
"""

import operator
import random
import sys
import warnings
import pandas as pd
import numpy as np
import scipy
from annoy import AnnoyIndex
from intervaltree import IntervalTree
from fbpca import pca
from scipy.sparse import csc_matrix, csr_matrix, vstack
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize

# Default parameters.
ALPHA = 0.10
APPROX = True
BATCH_SIZE = 5000
DIMRED = 100
HVG = None
KNN = 20
N_ITER = 500
PERPLEXITY = 1200
SIGMA = 15
VERBOSE = 2


def batch_correct_pc(adata, batch_series, correct=False,
                     n_components=30, sigma=25, alpha=0.1, knn=30, metric='angular',
                     **scanorama_kws):
    """
    Batch correction PCA based on integration

    Parameters
    ----------
    adata
        one major adata
    batch_series
        batch_series used for splitting adata
    correct
        if True, adata.X will be corrected inplace, otherwise only corrected PCs are added to adata.obsm['X_pca']
    n_components
        number of components in PCA
    sigma
        Correction smoothing parameter on Gaussian kernel.
    alpha
        Alignment score minimum cutoff.
    knn
        Number of nearest neighbors to use for matching.
    metric
        Metric to use in calculating KNN
    scanorama_kws
        Other Parameters passed to integration function
    Returns
    -------
    adata
    """

    scanorama_kws['dimred'] = n_components
    scanorama_kws['sigma'] = sigma
    scanorama_kws['alpha'] = alpha
    scanorama_kws['knn'] = knn
    scanorama_kws['metric'] = metric

    adata.obs['batch'] = batch_series
    adata_list = []
    indexes = []
    for _, sub_df in adata.obs.groupby('batch'):
        adata_list.append(adata[sub_df.index, :])
        indexes.extend(sub_df.index.tolist())

    if correct:
        integrated, corrected = correct_scanpy(adata_list, return_dimred=True,
                                                         **scanorama_kws)
        adata.X = np.vstack([ann.X.toarray() for ann in corrected])
    else:
        integrated = integrate_scanpy(adata_list, **scanorama_kws)

    pca_df = pd.DataFrame(np.vstack(integrated), index=indexes).reindex(adata.obs_names)

    adata.obsm['X_pca'] = pca_df.values
    # TODO fill up other PC related parts same as sc.tl.pca
    return adata

# Do batch correction on a list of data sets.
def correct(datasets_full, genes_list, return_dimred=False,
            batch_size=BATCH_SIZE, verbose=VERBOSE, ds_names=None,
            dimred=DIMRED, approx=APPROX, sigma=SIGMA, alpha=ALPHA, knn=KNN,
            return_dense=False, hvg=None, union=False,
            geosketch=False, geosketch_max=20000, seed=0, metric='manhattan'):
    """Integrate and batch correct a list of data sets.

    Parameters
    ----------
    datasets_full : `list` of `scipy.sparse.csr_matrix` or of `numpy.ndarray`
        Data sets to integrate and correct.
    genes_list: `list` of `list` of `string`
        List of genes for each data set.
    return_dimred: `bool`, optional (default: `False`)
        In addition to returning batch corrected matrices, also returns
        integrated low-dimesional embeddings.
    batch_size: `int`, optional (default: `5000`)
        The batch size used in the alignment vector computation. Useful when
        correcting very large (>100k samples) data sets. Set to large value
        that runs within available memory.
    verbose: `bool` or `int`, optional (default: 2)
        When `True` or not equal to 0, prints logging output.
    ds_names: `list` of `string`, optional
        When `verbose=True`, reports data set names in logging output.
    dimred: `int`, optional (default: 100)
        Dimensionality of integrated embedding.
    approx: `bool`, optional (default: `True`)
        Use approximate nearest neighbors, greatly speeds up matching runtime.
    sigma: `float`, optional (default: 15)
        Correction smoothing parameter on Gaussian kernel.
    alpha: `float`, optional (default: 0.10)
        Alignment score minimum cutoff.
    knn: `int`, optional (default: 20)
        Number of nearest neighbors to use for matching.
    return_dense: `bool`, optional (default: `False`)
        Return `numpy.ndarray` matrices instead of `scipy.sparse.csr_matrix`.
    hvg: `int`, optional (default: None)
        Use this number of top highly variable genes based on dispersion.
    seed: `int`, optional (default: 0)
        Random seed to use.

    Returns
    -------
    corrected, genes
        By default (`return_dimred=False`), returns a two-tuple containing a
        list of `scipy.sparse.csr_matrix` each with batch corrected values,
        and a single list of genes containing the intersection of inputted
        genes.

    integrated, corrected, genes
        When `return_dimred=False`, returns a three-tuple containing a list
        of `numpy.ndarray` with integrated low dimensional embeddings, a list
        of `scipy.sparse.csr_matrix` each with batch corrected values, and a
        a single list of genes containing the intersection of inputted genes.
    """
    np.random.seed(seed)
    random.seed(seed)

    datasets_full = check_datasets(datasets_full)

    datasets, genes = merge_datasets(datasets_full, genes_list,
                                     ds_names=ds_names, union=union)
    datasets_dimred, genes = process_data(datasets, genes, hvg=hvg,
                                          dimred=dimred)

    datasets_dimred = assemble(
        datasets_dimred,  # Assemble in low dimensional space.
        expr_datasets=datasets,  # Modified in place.
        verbose=verbose, knn=knn, sigma=sigma, approx=approx,
        alpha=alpha, ds_names=ds_names, batch_size=batch_size,
        geosketch=geosketch, geosketch_max=geosketch_max, metric=metric
    )

    if return_dense:
        datasets = [ds.toarray() for ds in datasets]

    if return_dimred:
        return datasets_dimred, datasets, genes

    return datasets, genes


# Integrate a list of data sets.
def integrate(datasets_full, genes_list, batch_size=BATCH_SIZE,
              verbose=VERBOSE, ds_names=None, dimred=DIMRED, approx=APPROX,
              sigma=SIGMA, alpha=ALPHA, knn=KNN, geosketch=False,
              geosketch_max=20000, n_iter=1, union=False, hvg=None, seed=0,
              metric='manhattan'):
    """Integrate a list of data sets.

    Parameters
    ----------
    datasets_full : `list` of `scipy.sparse.csr_matrix` or of `numpy.ndarray`
        Data sets to integrate and correct.
    genes_list: `list` of `list` of `string`
        List of genes for each data set.
    batch_size: `int`, optional (default: `5000`)
        The batch size used in the alignment vector computation. Useful when
        correcting very large (>100k samples) data sets. Set to large value
        that runs within available memory.
    verbose: `bool` or `int`, optional (default: 2)
        When `True` or not equal to 0, prints logging output.
    ds_names: `list` of `string`, optional
        When `verbose=True`, reports data set names in logging output.
    dimred: `int`, optional (default: 100)
        Dimensionality of integrated embedding.
    approx: `bool`, optional (default: `True`)
        Use approximate nearest neighbors, greatly speeds up matching runtime.
    sigma: `float`, optional (default: 15)
        Correction smoothing parameter on Gaussian kernel.
    alpha: `float`, optional (default: 0.10)
        Alignment score minimum cutoff.
    knn: `int`, optional (default: 20)
        Number of nearest neighbors to use for matching.
    hvg: `int`, optional (default: None)
        Use this number of top highly variable genes based on dispersion.
    seed: `int`, optional (default: 0)
        Random seed to use.

    Returns
    -------
    integrated, genes
        Returns a two-tuple containing a list of `numpy.ndarray` with
        integrated low dimensional embeddings and a single list of genes
        containing the intersection of inputted genes.
    """
    np.random.seed(seed)
    random.seed(seed)

    datasets_full = check_datasets(datasets_full)

    datasets, genes = merge_datasets(datasets_full, genes_list,
                                     ds_names=ds_names, union=union)
    datasets_dimred, genes = process_data(datasets, genes, hvg=hvg,
                                          dimred=dimred)

    for _ in range(n_iter):
        datasets_dimred = assemble(
            datasets_dimred,  # Assemble in low dimensional space.
            verbose=verbose, knn=knn, sigma=sigma, approx=approx,
            alpha=alpha, ds_names=ds_names, batch_size=batch_size,
            geosketch=geosketch, geosketch_max=geosketch_max, metric=metric,
        )

    return datasets_dimred, genes


# Batch correction with scanpy's AnnData object.
def correct_scanpy(adatas, **kwargs):
    """Batch correct a list of `scanpy.api.AnnData`.

    Parameters
    ----------
    adatas : `list` of `scanpy.api.AnnData`
        Data sets to integrate and/or correct.
    kwargs : `dict`
        See documentation for the `correct()` method for a full list of
        parameters to use for batch correction.

    Returns
    -------
    corrected
        By default (`return_dimred=False`), returns a list of new
        `scanpy.api.AnnData`.

    integrated, corrected
        When `return_dimred=True`, returns a two-tuple containing a list of
        `np.ndarray` with integrated low-dimensional embeddings and a list
        of new `scanpy.api.AnnData`.
    """
    if 'return_dimred' in kwargs and kwargs['return_dimred']:
        datasets_dimred, datasets, genes = correct(
            [adata.X for adata in adatas],
            [adata.var_names.values for adata in adatas],
            **kwargs
        )
    else:
        datasets, genes = correct(
            [adata.X for adata in adatas],
            [adata.var_names.values for adata in adatas],
            **kwargs
        )

    from anndata import AnnData

    new_adatas = []
    for i in range(len((adatas))):
        adata = AnnData(datasets[i])
        adata.var_names = genes
        new_adatas.append(adata)

    if 'return_dimred' in kwargs and kwargs['return_dimred']:
        return datasets_dimred, new_adatas
    else:
        return new_adatas


# Integration with scanpy's AnnData object.
def integrate_scanpy(adatas, **kwargs):
    """Integrate a list of `scanpy.api.AnnData`.

    Parameters
    ----------
    adatas : `list` of `scanpy.api.AnnData`
        Data sets to integrate.
    kwargs : `dict`
        See documentation for the `integrate()` method for a full list of
        parameters to use for batch correction.

    Returns
    -------
    integrated
        Returns a list of `np.ndarray` with integrated low-dimensional
        embeddings.
    """
    datasets_dimred, genes = integrate(
        [adata.X for adata in adatas],
        [adata.var_names.values for adata in adatas],
        **kwargs
    )

    return datasets_dimred


# Put datasets into a single matrix with the intersection of all genes.
def merge_datasets(datasets, genes, ds_names=None, verbose=True,
                   union=False):
    if union:
        sys.stderr.write(
            'WARNING: Integrating based on the union of genes is '
            'highly discouraged, consider taking the intersection '
            'or requantifying gene expression.\n'
        )

    # Find genes in common.
    keep_genes = set()
    for idx, gene_list in enumerate(genes):
        if len(keep_genes) == 0:
            keep_genes = set(gene_list)
        elif union:
            keep_genes |= set(gene_list)
        else:
            keep_genes &= set(gene_list)
        if not union and not ds_names is None and verbose:
            print('After {}: {} genes'.format(ds_names[idx], len(keep_genes)))
        if len(keep_genes) == 0:
            print('Error: No genes found in all datasets, exiting...')
            exit(1)
    if verbose:
        print('Found {} genes among all datasets'
              .format(len(keep_genes)))

    if union:
        union_genes = sorted(keep_genes)
        for i in range(len(datasets)):
            if verbose:
                print('Processing data set {}'.format(i))
            X_new = np.zeros((datasets[i].shape[0], len(union_genes)))
            X_old = csc_matrix(datasets[i])
            gene_to_idx = {gene: idx for idx, gene in enumerate(genes[i])}
            for j, gene in enumerate(union_genes):
                if gene in gene_to_idx:
                    X_new[:, j] = X_old[:, gene_to_idx[gene]].toarray().flatten()
            datasets[i] = csr_matrix(X_new)
        ret_genes = np.array(union_genes)
    else:
        # Only keep genes in common.
        ret_genes = np.array(sorted(keep_genes))
        for i in range(len(datasets)):
            # Remove duplicate genes.
            uniq_genes, uniq_idx = np.unique(genes[i], return_index=True)
            datasets[i] = datasets[i][:, uniq_idx]

            # Do gene filtering.
            gene_sort_idx = np.argsort(uniq_genes)
            gene_idx = [idx for idx in gene_sort_idx
                        if uniq_genes[idx] in keep_genes]
            datasets[i] = datasets[i][:, gene_idx]
            assert (np.array_equal(uniq_genes[gene_idx], ret_genes))

    return datasets, ret_genes


def check_datasets(datasets_full):
    datasets_new = []
    for i, ds in enumerate(datasets_full):
        if issubclass(type(ds), np.ndarray):
            datasets_new.append(csr_matrix(ds))
        elif issubclass(type(ds), scipy.sparse.csr.csr_matrix):
            datasets_new.append(ds)
        else:
            sys.stderr.write('ERROR: Data sets must be numpy array or '
                             'scipy.sparse.csr_matrix, received type '
                             '{}.\n'.format(type(ds)))
            exit(1)
    return datasets_new


def reduce_dimensionality(X, dim_red_k=100):
    k = min((dim_red_k, X.shape[0], X.shape[1]))
    U, s, Vt = pca(X, k=k)  # Automatically centers.
    return U[:, range(k)] * s[range(k)]


# Randomized SVD.
def dimensionality_reduce(datasets, dimred=DIMRED):
    X = vstack(datasets)
    X = reduce_dimensionality(X, dim_red_k=dimred)
    datasets_dimred = []
    base = 0
    for ds in datasets:
        datasets_dimred.append(X[base:(base + ds.shape[0]), :])
        base += ds.shape[0]
    return datasets_dimred


def dispersion(X):
    mean = X.mean(0)
    dispersion = np.zeros(mean.shape)
    nonzero_idx = np.nonzero(mean > 1e-10)[1]
    X_nonzero = X[:, nonzero_idx]
    nonzero_mean = X_nonzero.mean(0)
    nonzero_var = (X_nonzero.multiply(X_nonzero)).mean(0)
    temp = (nonzero_var / nonzero_mean)
    dispersion[mean > 1e-10] = temp.A1
    dispersion[mean <= 1e-10] = float('-inf')
    return dispersion


# Normalize and reduce dimensionality.
def process_data(datasets, genes, hvg=HVG, dimred=DIMRED, verbose=True):
    # Only keep highly variable genes
    if not hvg is None and hvg > 0 and hvg < len(genes):
        if verbose:
            print('Highly variable filter...')
        X = vstack(datasets)
        disp = dispersion(X)
        highest_disp_idx = np.argsort(disp[0])[::-1]
        top_genes = set(genes[highest_disp_idx[range(hvg)]])
        for i in range(len(datasets)):
            gene_idx = [idx for idx, g_i in enumerate(genes)
                        if g_i in top_genes]
            datasets[i] = datasets[i][:, gene_idx]
        genes = np.array(sorted(top_genes))

    # Normalize.
    if verbose:
        print('Normalizing...')
    for i, ds in enumerate(datasets):
        datasets[i] = normalize(ds, axis=1)

    # Compute compressed embedding.
    if dimred > 0:
        if verbose:
            print('Reducing dimension...')
        datasets_dimred = dimensionality_reduce(datasets, dimred=dimred)
        if verbose:
            print('Done processing.')
        return datasets_dimred, genes

    if verbose:
        print('Done processing.')

    return datasets, genes


# Exact nearest neighbors search.
def nn(ds1, ds2, knn=KNN, metric_p=2):
    # Find nearest neighbors of first dataset.
    nn_ = NearestNeighbors(knn, p=metric_p)
    nn_.fit(ds2)
    ind = nn_.kneighbors(ds1, return_distance=False)

    match = set()
    for a, b in zip(range(ds1.shape[0]), ind):
        for b_i in b:
            match.add((a, b_i))

    return match


# Approximate nearest neighbors using locality sensitive hashing.
def nn_approx(ds1, ds2, knn=KNN, metric='manhattan', n_trees=10):
    # Build index.
    print(f'AnnoyIndex using {metric} metric')
    a = AnnoyIndex(ds2.shape[1], metric=metric)
    for i in range(ds2.shape[0]):
        a.add_item(i, ds2[i, :])
    a.build(n_trees)

    # Search index.
    ind = []
    for i in range(ds1.shape[0]):
        ind.append(a.get_nns_by_vector(ds1[i, :], knn, search_k=-1))
    ind = np.array(ind)

    # Match.
    match = set()
    for a, b in zip(range(ds1.shape[0]), ind):
        for b_i in b:
            match.add((a, b_i))

    return match


# Populate a table (in place) that stores mutual nearest neighbors
# between datasets.
def fill_table(table, i, curr_ds, datasets, base_ds=0,
               knn=KNN, approx=APPROX, metric='manhattan'):
    curr_ref = np.concatenate(datasets)
    if approx:
        match = nn_approx(curr_ds, curr_ref, knn=knn, metric=metric)
    else:
        match = nn(curr_ds, curr_ref, knn=knn, metric_p=1)

    # Build interval tree.
    itree_ds_idx = IntervalTree()
    itree_pos_base = IntervalTree()
    pos = 0
    for j in range(len(datasets)):
        n_cells = datasets[j].shape[0]
        itree_ds_idx[pos:(pos + n_cells)] = base_ds + j
        itree_pos_base[pos:(pos + n_cells)] = pos
        pos += n_cells

    # Store all mutual nearest neighbors between datasets.
    for d, r in match:
        interval = itree_ds_idx[r]
        assert (len(interval) == 1)
        j = interval.pop().data
        interval = itree_pos_base[r]
        assert (len(interval) == 1)
        base = interval.pop().data
        if not (i, j) in table:
            table[(i, j)] = set()
        table[(i, j)].add((d, r - base))
        assert (r - base >= 0)


gs_idxs = None


# Fill table of alignment scores.
def find_alignments_table(datasets, knn=KNN, approx=APPROX, verbose=VERBOSE,
                          prenormalized=False, geosketch=False,
                          geosketch_max=20000, metric='manhattan'):
    if not prenormalized:
        datasets = [normalize(ds, axis=1) for ds in datasets]

    if geosketch:
        # Only match cells in geometric sketches.
        raise NotImplementedError

    table = {}
    for i in range(len(datasets)):
        if len(datasets[:i]) > 0:
            fill_table(table, i, datasets[i], datasets[:i], knn=knn,
                       approx=approx, metric=metric)
        if len(datasets[i + 1:]) > 0:
            fill_table(table, i, datasets[i], datasets[i + 1:],
                       knn=knn, base_ds=i + 1, approx=approx, metric=metric)
    # Count all mutual nearest neighbors between datasets.
    matches = {}
    table1 = {}
    if verbose > 1:
        table_print = np.zeros((len(datasets), len(datasets)))
    for i in range(len(datasets)):
        for j in range(len(datasets)):
            if i >= j:
                continue
            if not (i, j) in table or not (j, i) in table:
                continue
            match_ij = table[(i, j)]
            match_ji = set([(b, a) for a, b in table[(j, i)]])
            matches[(i, j)] = match_ij & match_ji

            table1[(i, j)] = (max(
                float(len(set([idx for idx, _ in matches[(i, j)]]))) /
                datasets[i].shape[0],
                float(len(set([idx for _, idx in matches[(i, j)]]))) /
                datasets[j].shape[0]
            ))
            if verbose > 1:
                table_print[i, j] += table1[(i, j)]

            if geosketch:
                # Translate matches within geometric sketches to original indices.
                matches_mnn = matches[(i, j)]
                matches[(i, j)] = [
                    (gs_idxs[i][a], gs_idxs[j][b]) for a, b in matches_mnn
                ]

    if verbose > 1:
        print(table_print)
        return table1, table_print, matches
    else:
        return table1, None, matches


# Find the matching pairs of cells between datasets.
def find_alignments(datasets, knn=KNN, approx=APPROX, verbose=VERBOSE,
                    alpha=ALPHA, prenormalized=False,
                    geosketch=False, geosketch_max=20000, metric='manhattan'):
    table1, _, matches = find_alignments_table(
        datasets, knn=knn, approx=approx, verbose=verbose,
        prenormalized=prenormalized,
        geosketch=geosketch, geosketch_max=geosketch_max,
        metric=metric
    )

    alignments = [(i, j) for (i, j), val in reversed(
        sorted(table1.items(), key=operator.itemgetter(1))
    ) if val > alpha]

    return alignments, matches


# Find connections between datasets to identify panoramas.
def connect(datasets, knn=KNN, approx=APPROX, alpha=ALPHA,
            verbose=VERBOSE, metric='manhattan'):
    # Find alignments.
    alignments, _ = find_alignments(
        datasets, knn=knn, approx=approx, alpha=alpha,
        verbose=verbose, metric=metric,
    )
    if verbose:
        print(alignments)

    panoramas = []
    connected = set()
    for i, j in alignments:
        # See if datasets are involved in any current panoramas.
        panoramas_i = [panoramas[p] for p in range(len(panoramas))
                       if i in panoramas[p]]
        assert (len(panoramas_i) <= 1)
        panoramas_j = [panoramas[p] for p in range(len(panoramas))
                       if j in panoramas[p]]
        assert (len(panoramas_j) <= 1)

        if len(panoramas_i) == 0 and len(panoramas_j) == 0:
            panoramas.append([i])
            panoramas_i = [panoramas[-1]]

        if len(panoramas_i) == 0:
            panoramas_j[0].append(i)
        elif len(panoramas_j) == 0:
            panoramas_i[0].append(j)
        elif panoramas_i[0] != panoramas_j[0]:
            panoramas_i[0] += panoramas_j[0]
            panoramas.remove(panoramas_j[0])

        connected.add(i)
        connected.add(j)

    for i in range(len(datasets)):
        if not i in connected:
            panoramas.append([i])

    return panoramas


def handle_zeros_in_scale(scale, copy=True):
    ''' Makes sure that whenever scale is zero, we handle it correctly.
    This happens in most scalers when we have constant features.
    Adapted from sklearn.preprocessing.data'''

    # if we are fitting on 1D arrays, scale might be a scalar
    if np.isscalar(scale):
        if scale == .0:
            scale = 1.
        return scale
    elif isinstance(scale, np.ndarray):
        if copy:
            # New array to avoid side-effects
            scale = scale.copy()
        scale[scale == 0.0] = 1.0
    return scale


# To reduce memory usage, split bias computation into batches.
def batch_bias(curr_ds, match_ds, bias, batch_size=None, sigma=SIGMA):
    if batch_size is None:
        weights = rbf_kernel(curr_ds, match_ds, gamma=0.5 * sigma)
        weights = normalize(weights, axis=1, norm='l1')
        avg_bias = np.dot(weights, bias)
        return avg_bias

    base = 0
    avg_bias = np.zeros(curr_ds.shape)
    denom = np.zeros(curr_ds.shape[0])
    while base < match_ds.shape[0]:
        batch_idx = range(
            base, min(base + batch_size, match_ds.shape[0])
        )
        weights = rbf_kernel(curr_ds, match_ds[batch_idx, :],
                             gamma=0.5 * sigma)
        avg_bias += np.dot(weights, bias[batch_idx, :])
        denom += np.sum(weights, axis=1)
        base += batch_size

    denom = handle_zeros_in_scale(denom, copy=False)
    avg_bias /= denom[:, np.newaxis]

    return avg_bias


# Compute nonlinear translation vectors between dataset
# and a reference.
def transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma=SIGMA, cn=False,
              batch_size=None):
    # Compute the matching.
    match_ds = curr_ds[ds_ind, :]
    match_ref = curr_ref[ref_ind, :]
    bias = match_ref - match_ds
    if cn:
        match_ds = match_ds.toarray()
        curr_ds = curr_ds.toarray()
        bias = bias.toarray()

    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            avg_bias = batch_bias(curr_ds, match_ds, bias, sigma=sigma,
                                  batch_size=batch_size)
        except RuntimeWarning:
            sys.stderr.write('WARNING: Oversmoothing detected, refusing to batch '
                             'correct, consider lowering sigma value.\n')
            return csr_matrix(curr_ds.shape, dtype=float)
        except MemoryError:
            if batch_size is None:
                sys.stderr.write('WARNING: Out of memory, consider turning on '
                                 'batched computation with batch_size parameter.\n')
            else:
                sys.stderr.write('WARNING: Out of memory, consider lowering '
                                 'the batch_size parameter.\n')
            return csr_matrix(curr_ds.shape, dtype=float)

    if cn:
        avg_bias = csr_matrix(avg_bias)

    return avg_bias


# Finds alignments between datasets and uses them to construct
# panoramas. "Merges" datasets by correcting gene expression
# values.
def assemble(datasets, verbose=VERBOSE, knn=KNN,
             sigma=SIGMA, approx=APPROX, alpha=ALPHA, expr_datasets=None,
             ds_names=None, batch_size=None,
             geosketch=False, geosketch_max=20000, alignments=None, matches=None,
             metric='manhattan'):
    if len(datasets) == 1:
        return datasets

    if alignments is None and matches is None:
        alignments, matches = find_alignments(
            datasets, knn=knn, approx=approx, alpha=alpha, verbose=verbose,
            geosketch=geosketch, geosketch_max=geosketch_max, metric=metric
        )

    ds_assembled = {}
    panoramas = []
    for i, j in alignments:
        if verbose:
            if ds_names is None:
                print('Processing datasets {}'.format((i, j)))
            else:
                print('Processing datasets {} <=> {}'.
                      format(ds_names[i], ds_names[j]))

        # Only consider a dataset a fixed amount of times.
        if not i in ds_assembled:
            ds_assembled[i] = 0
        ds_assembled[i] += 1
        if not j in ds_assembled:
            ds_assembled[j] = 0
        ds_assembled[j] += 1
        if ds_assembled[i] > 3 and ds_assembled[j] > 3:
            continue

        # See if datasets are involved in any current panoramas.
        panoramas_i = [panoramas[p] for p in range(len(panoramas))
                       if i in panoramas[p]]
        assert (len(panoramas_i) <= 1)
        panoramas_j = [panoramas[p] for p in range(len(panoramas))
                       if j in panoramas[p]]
        assert (len(panoramas_j) <= 1)

        if len(panoramas_i) == 0 and len(panoramas_j) == 0:
            if datasets[i].shape[0] < datasets[j].shape[0]:
                i, j = j, i
            panoramas.append([i])
            panoramas_i = [panoramas[-1]]

        # Map dataset i to panorama j.
        if len(panoramas_i) == 0:
            curr_ds = datasets[i]
            curr_ref = np.concatenate([datasets[p] for p in panoramas_j[0]])

            match = []
            base = 0
            for p in panoramas_j[0]:
                if i < p and (i, p) in matches:
                    match.extend([(a, b + base) for a, b in matches[(i, p)]])
                elif i > p and (p, i) in matches:
                    match.extend([(b, a + base) for a, b in matches[(p, i)]])
                base += datasets[p].shape[0]

            ds_ind = [a for a, _ in match]
            ref_ind = [b for _, b in match]

            bias = transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma=sigma,
                             batch_size=batch_size)
            datasets[i] = curr_ds + bias

            if expr_datasets:
                curr_ds = expr_datasets[i]
                curr_ref = vstack([expr_datasets[p]
                                   for p in panoramas_j[0]])
                bias = transform(curr_ds, curr_ref, ds_ind, ref_ind,
                                 sigma=sigma, cn=True, batch_size=batch_size)
                expr_datasets[i] = curr_ds + bias

            panoramas_j[0].append(i)

        # Map dataset j to panorama i.
        elif len(panoramas_j) == 0:
            curr_ds = datasets[j]
            curr_ref = np.concatenate([datasets[p] for p in panoramas_i[0]])

            match = []
            base = 0
            for p in panoramas_i[0]:
                if j < p and (j, p) in matches:
                    match.extend([(a, b + base) for a, b in matches[(j, p)]])
                elif j > p and (p, j) in matches:
                    match.extend([(b, a + base) for a, b in matches[(p, j)]])
                base += datasets[p].shape[0]

            ds_ind = [a for a, _ in match]
            ref_ind = [b for _, b in match]

            bias = transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma=sigma,
                             batch_size=batch_size)
            datasets[j] = curr_ds + bias

            if expr_datasets:
                curr_ds = expr_datasets[j]
                curr_ref = vstack([expr_datasets[p]
                                   for p in panoramas_i[0]])
                bias = transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma=sigma,
                                 cn=True, batch_size=batch_size)
                expr_datasets[j] = curr_ds + bias

            panoramas_i[0].append(j)

        # Merge two panoramas together.
        else:
            curr_ds = np.concatenate([datasets[p] for p in panoramas_i[0]])
            curr_ref = np.concatenate([datasets[p] for p in panoramas_j[0]])

            # Find base indices into each panorama.
            base_i = 0
            for p in panoramas_i[0]:
                if p == i: break
                base_i += datasets[p].shape[0]
            base_j = 0
            for p in panoramas_j[0]:
                if p == j: break
                base_j += datasets[p].shape[0]

            # Find matching indices.
            match = []
            base = 0
            for p in panoramas_i[0]:
                if p == i and j < p and (j, p) in matches:
                    match.extend([(b + base, a + base_j)
                                  for a, b in matches[(j, p)]])
                elif p == i and j > p and (p, j) in matches:
                    match.extend([(a + base, b + base_j)
                                  for a, b in matches[(p, j)]])
                base += datasets[p].shape[0]
            base = 0
            for p in panoramas_j[0]:
                if p == j and i < p and (i, p) in matches:
                    match.extend([(a + base_i, b + base)
                                  for a, b in matches[(i, p)]])
                elif p == j and i > p and (p, i) in matches:
                    match.extend([(b + base_i, a + base)
                                  for a, b in matches[(p, i)]])
                base += datasets[p].shape[0]

            ds_ind = [a for a, _ in match]
            ref_ind = [b for _, b in match]

            # Apply transformation to entire panorama.
            bias = transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma=sigma,
                             batch_size=batch_size)
            curr_ds += bias
            base = 0
            for p in panoramas_i[0]:
                n_cells = datasets[p].shape[0]
                datasets[p] = curr_ds[base:(base + n_cells), :]
                base += n_cells

            if not expr_datasets is None:
                curr_ds = vstack([expr_datasets[p]
                                  for p in panoramas_i[0]])
                curr_ref = vstack([expr_datasets[p]
                                   for p in panoramas_j[0]])
                bias = transform(curr_ds, curr_ref, ds_ind, ref_ind,
                                 sigma=sigma, cn=True, batch_size=batch_size)
                curr_ds += bias
                base = 0
                for p in panoramas_i[0]:
                    n_cells = expr_datasets[p].shape[0]
                    expr_datasets[p] = curr_ds[base:(base + n_cells), :]
                    base += n_cells

            # Merge panoramas i and j and delete one.
            if panoramas_i[0] != panoramas_j[0]:
                panoramas_i[0] += panoramas_j[0]
                panoramas.remove(panoramas_j[0])
    return datasets


def interpret_alignments(datasets, expr_datasets, genes,
                         verbose=VERBOSE, knn=KNN, approx=APPROX,
                         alpha=ALPHA, n_permutations=None, metric='manhattan'):
    if n_permutations is None:
        n_permutations = float(len(genes) * 30)

    alignments, matches = find_alignments(
        datasets, knn=knn, approx=approx, alpha=alpha, verbose=verbose, metric=metric
    )

    for i, j in alignments:
        # Compute average bias vector that aligns two datasets together.
        ds_i = expr_datasets[i]
        ds_j = expr_datasets[j]
        if i < j:
            match = matches[(i, j)]
        else:
            match = matches[(j, i)]
        i_ind = [a for a, _ in match]
        j_ind = [b for _, b in match]
        avg_bias = np.absolute(
            np.mean(ds_j[j_ind, :] - ds_i[i_ind, :], axis=0)
        )

        # Construct null distribution and compute p-value.
        null_bias = (
                ds_j[np.random.randint(ds_j.shape[0], size=n_permutations), :] -
                ds_i[np.random.randint(ds_i.shape[0], size=n_permutations), :]
        )
        p = ((np.sum(np.greater_equal(
            np.absolute(np.tile(avg_bias, (n_permutations, 1))),
            np.absolute(null_bias)
        ), axis=0, dtype=float) + 1) / (n_permutations + 1))

        print('>>>> Stats for alignment {}'.format((i, j)))
        for k in range(len(p)):
            print('{}\t{}'.format(genes[k], p[k]))

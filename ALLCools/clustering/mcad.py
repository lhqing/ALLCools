import numpy as np
from sklearn.decomposition import TruncatedSVD
import pandas as pd
import warnings
from pybedtools import BedTool


def remove_black_list_region(adata, black_list_path, f=0.2):
    """
    Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path

    Parameters
    ----------
    adata
    black_list_path
        Path to the black list bed file
    f
        Fraction of overlap when calling bedtools intersect
    Returns
    -------
    None
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_bed_df = adata.var[['chrom', 'start', 'end']]
        feature_bed = BedTool.from_dataframe(feature_bed_df)
        black_list_bed = BedTool(black_list_path)
        black_feature = feature_bed.intersect(black_list_bed, f=f, wa=True)
        black_feature_index = black_feature.to_dataframe().set_index(
            ['chrom', 'start', 'end']).index
    black_feature_id = pd.Index(
        feature_bed_df.reset_index().set_index(['chrom', 'start', 'end']).loc[black_feature_index]['region'])

    print(f'{black_feature_id.size} features removed due to overlapping'
          f' (bedtools intersect -f {f}) with black list regions.')
    adata._inplace_subset_var(~adata.var_names.isin(black_feature_id))
    return


def binarize_matrix(adata, cutoff=0.95):
    """
    Binarize adata.X with adata.X > cutoff

    Parameters
    ----------
    adata
        AnnData object whose X is survival function matrix
    cutoff
        Cutoff to binarize the survival function

    Returns
    -------
    None
    """
    adata.X = (adata.X > cutoff).astype(np.int8)
    return


def filter_regions(adata, hypo_cutoff=None):
    """
    Filter regions based on their

    Parameters
    ----------
    adata
    hypo_cutoff
        min number of cells that are hypo-methylated (1) in this region.
        If None, will use adata.shape[0] * 0.003

    Returns
    -------

    """
    if hypo_cutoff is None:
        hypo_cutoff = adata.shape[0] * 0.003
    hypo_judge = adata.X.sum(axis=0).A1 > hypo_cutoff
    adata._inplace_subset_var(hypo_judge)
    return


def lsi(adata, scale_factor=100000, n_components=50, algorithm='arpack', obsm='X_pca', random_state=0):
    """
    Run TF-IDF on the binarized adata.X, followed by TruncatedSVD and then scale the components by svd.singular_values_

    Parameters
    ----------
    adata
    scale_factor
    n_components
    algorithm
    obsm
    random_state

    Returns
    -------

    """
    if n_components == 'auto':
        n_components = 100
    # tf-idf
    data = adata.X.astype(np.int8).copy()
    col_sum = data.sum(axis=0).A1
    row_sum = data.sum(axis=1).A1

    idf = np.log(1 + data.shape[0] / col_sum)
    tf = data
    tf.data = tf.data / np.repeat(row_sum, row_sum)
    tf.data = np.log(tf.data * scale_factor + 1)
    tf = tf.multiply(idf)
    svd = TruncatedSVD(n_components=n_components, algorithm=algorithm, random_state=random_state)
    matrix_reduce = svd.fit_transform(tf)
    matrix_reduce = matrix_reduce / svd.singular_values_

    # PCA is the default name for many following steps in scanpy, use the name here for convenience.
    # However, this is not PCA
    adata.obsm[obsm] = matrix_reduce
    return

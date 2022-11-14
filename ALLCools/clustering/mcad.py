import warnings

import numpy as np
import pandas as pd
from pybedtools import BedTool


def remove_black_list_region(adata, black_list_path, region_axis=1, f=0.2):
    """
    Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path.

    Parameters
    ----------
    adata
        AnnData object
    black_list_path
        Path to the black list bed file
    region_axis
        Axis of regions. 0 for adata.obs, 1 for adata.var
    f
        Fraction of overlap when calling bedtools intersect
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if region_axis == 1:
            feature_bed_df = adata.var[["chrom", "start", "end"]]
        elif region_axis == 0:
            feature_bed_df = adata.obs[["chrom", "start", "end"]]
        else:
            raise ValueError("region_axis should be 0 or 1.")
        feature_bed = BedTool.from_dataframe(feature_bed_df)
        black_list_bed = BedTool(black_list_path)
        black_feature = feature_bed.intersect(black_list_bed, f=f, wa=True)
        try:
            black_feature_index = black_feature.to_dataframe().set_index(["chrom", "start", "end"]).index
            black_feature_id = pd.Index(
                feature_bed_df.reset_index()
                .set_index(["chrom", "start", "end"])
                .loc[black_feature_index][feature_bed_df.index.name]
            )
            print(
                f"{black_feature_id.size} regions removed due to overlapping"
                f" (bedtools intersect -f {f}) with black list regions."
            )
            if region_axis == 1:
                adata._inplace_subset_var(~adata.var_names.isin(black_feature_id))
            else:
                adata._inplace_subset_obs(~adata.obs_names.isin(black_feature_id))
        except pd.errors.EmptyDataError:
            # no overlap with black list
            pass
    return


def remove_chromosomes(adata, exclude_chromosomes=None, include_chromosomes=None, chrom_col="chrom"):
    """Remove chromosomes from adata.var."""
    judge = None
    if exclude_chromosomes is not None:
        not_to_exclude = ~adata.var[chrom_col].isin(exclude_chromosomes)
        judge = not_to_exclude
    if include_chromosomes is not None:
        include = adata.var[chrom_col].isin(include_chromosomes)
        if judge is None:
            judge = include
        else:
            judge &= include

    if judge is not None:
        adata._inplace_subset_var(judge)
        print(f"{adata.shape[1]} regions remained.")
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


def filter_regions(adata, hypo_percent=0.5, n_cell=None, zscore_abs_cutoff=None):
    """
    Filter regions based on % of cells having non-zero scores.

    Parameters
    ----------
    adata
        AnnData object
    hypo_percent
        min % of cells that are non-zero in this region. If n_cell is provided, this parameter will be ignored.
    n_cell
        number of cells that are non-zero in this region.
    zscore_abs_cutoff
        absolute feature non-zero cell count zscore cutoff to remove lowest and highest coverage features.
    """
    _nnz = (adata.X > 0).sum(axis=0)
    try:
        feature_nnz_cell = _nnz.A1
    except AttributeError:
        feature_nnz_cell = _nnz.ravel()

    if n_cell is None:
        n_cell = int(adata.shape[0] * hypo_percent / 100)
    n_cell_judge = feature_nnz_cell > n_cell
    adata._inplace_subset_var(n_cell_judge)
    feature_nnz_cell = feature_nnz_cell[n_cell_judge].copy()

    if zscore_abs_cutoff is not None:
        from scipy.stats import zscore

        zscore_judge = np.abs(zscore(np.log2(feature_nnz_cell))) < zscore_abs_cutoff
        adata._inplace_subset_var(zscore_judge)

    print(f"{adata.shape[1]} regions remained.")
    return

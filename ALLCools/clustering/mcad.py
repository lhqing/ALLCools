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


def filter_regions(adata, hypo_percent=0.5, n_cell=None, zscore_abs_cutoff=None, chunk_size=5000):
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
    chunk_size: int
        chunk the regions (columns) using chunk_size.
    """

    _nnz = []
    ncols = adata.X.shape[1]
    for i in range(ncols // chunk_size + 1):
        _nnz.append(
            (adata.X[:, i * chunk_size : (i + 1) * chunk_size] > 0).sum(axis=0)
        )  # number of cells with non-hyper-methylation

    if n_cell is None:
        n_cell = int(adata.shape[0] * hypo_percent / 100)

    feature_nnz_cell = []
    for chunk in _nnz:
        try:
            a = chunk.A1
        except AttributeError:
            a = chunk.ravel()
        feature_nnz_cell.append(a)

    feature_nnz_cell = np.concatenate(feature_nnz_cell)  # in default, axis=0 for np.concatenate

    n_cell_judge = feature_nnz_cell > n_cell
    adata._inplace_subset_var(n_cell_judge)
    feature_nnz_cell = feature_nnz_cell[n_cell_judge].copy()

    if zscore_abs_cutoff is not None:
        from scipy.stats import zscore

        zscore_judge = np.abs(zscore(np.log2(feature_nnz_cell))) < zscore_abs_cutoff
        adata._inplace_subset_var(zscore_judge)

    print(f"{adata.shape[1]} regions remained.")
    return


def filter_regions_hyper(adata, non_hyper_percent=0.5, n_cell=None, zscore_abs_cutoff=None, chunk_size=5000):
    """
    Filter regions based on % of cells having non-zero scores.

    Parameters
    ----------
    adata
        AnnData object
    non_hyper_percent
        min % of cells that are zero in this region. If n_cell is provided, this parameter will be ignored.
    n_cell
        number of cells that are non-zero in this region.
    zscore_abs_cutoff
        absolute feature non-zero cell count zscore cutoff to remove lowest and highest coverage features.
    chunk_size: int
        chunk the columns (regions).
    """
    _nnz = []
    ncols = adata.X.shape[1]
    n_batch=ncols // chunk_size + 1
    print("Chunking regions..")
    for i in range(n_batch):
        print(round((i+1) / n_batch * 100, 2), '%\t\t\t', end='\r')
        _nnz.append(
            (adata.X[:, i * chunk_size : (i + 1) * chunk_size] == 0).sum(axis=0)
        )  # number of cells with non-hyper-methylation

    if n_cell is None:
        n_cell = int(adata.shape[0] * non_hyper_percent / 100)

    feature_nnz_cell = []
    for chunk in _nnz:
        try:
            a = chunk.A1
        except AttributeError:
            a = chunk.ravel()
        feature_nnz_cell.append(a)

    feature_nnz_cell = np.concatenate(feature_nnz_cell)  # in default, axis=0 for np.concatenate
    n_cell_judge = feature_nnz_cell > n_cell

    adata._inplace_subset_var(n_cell_judge)
    feature_nnz_cell = feature_nnz_cell[n_cell_judge].copy()

    if zscore_abs_cutoff is not None:
        from scipy.stats import zscore

        zscore_judge = np.abs(zscore(np.log2(feature_nnz_cell))) < zscore_abs_cutoff
        adata._inplace_subset_var(zscore_judge)

    print(f"{adata.shape[1]} regions remained.")
    return

def cal_chi2_pvalue(adata=None,df=None,cols=None,n_group_a=None,
                    n_group_b=None,comparison="GroupA | groupB",
                    p_cutoff=0.05,added_key='feature_chi2'):
    def get_pvalue(x):
        try:
            p = chi2_contingency(
            [
                [x['GroupA'],x['GroupB']],
                [x['NonGroupA'],x['NonGroupB']]
            ])[1]
        except:
            p=np.nan
        return p
    
    df_stat=df.groupby('Comparison')[cols].agg(np.sum).T
    df_stat['NonGroupA']=n_group_a-df_stat['GroupA']
    df_stat['NonGroupB']=n_group_b-df_stat['GroupB']
    df_stat['chi2_p']=df_stat.apply(lambda x:get_pvalue(x),axis=1)
    df_stat=df_stat.loc[~ df_stat.chi2_p.isna()]
    if df_stat.shape[0]==0:
        return None
    df_stat['Comparison']=comparison
    df_stat['chi2_adjp']=df_stat['chi2_p'] * df_stat.shape[0]
    df_stat=df_stat.loc[df_stat['chi2_adjp'] <= p_cutoff]
    if df_stat.shape[0]==0:
        return None
    if adata.uns[added_key] is None:
        adata.uns[added_key]=df_stat
    else:
        adata.uns[added_key]=pd.concat([adata.uns[added_key],df_stat])

def filter_region_chi2(adata,group=None,chunk_size=2000,
                         p_cutoff=0.05,added_key='feature_chi2',
                         pair_wise=True,topn=1500):
    """
    Filter regions based on Chi2 square test

    Parameters
    ----------
    adata: adata
    group: str
        group name to perform chi2 square test.
    chunk_size: int
        chunk columns (regions) using chunk_size
    p_cutoff: float
        p value cut off for chi2 square test adjusted pvalue.
    added_key: str
        key added to adata.uns
    pair_wise:
        In addition to one verse rest, whether to perform one verse one 
        comparison too.
    topn: int
        top n features to selected for each comparison. After top n features
        are selected for each comparison, the union features names will be the 
        final features, so the final number could be different than
        n_comparison * topn.
    """
    from scipy.stats import chi2_contingency
    import itertools
    assert not group is None
    ncols = adata.X.shape[1]
    nrows=adata.X.shape[0]
    n_batch=ncols // chunk_size + 1
    print("Chunking regions..")
    n_count_dict=adata.obs[group].value_counts().to_dict()
    adata.uns[added_key]=None
    for i in range(n_batch): #chunk columns (regions)
        print(round((i+1) / n_batch * 100, 2), '%\t\t\t', end='\r')
        cols=adata.var_names[i * chunk_size : (i + 1) * chunk_size]
        df=pd.DataFrame(adata.X[:, i * chunk_size : (i + 1) * chunk_size].toarray(),
                        columns=cols,index=adata.obs_names
        )
        df[group]=adata.obs[group]
        
        for g in df[group].unique():
            # one verse rest
            df['Comparison']=df[group].apply(lambda x:"GroupA" if x==g else 'GroupB')
            n_group_a=n_count_dict[g]
            n_group_b=nrows - n_group_a
            cal_chi2_pvalue(adata,df,cols=cols, n_group_a=n_group_a,
                            n_group_b=n_group_b,comparison=g+" | "+"Rest",
                            p_cutoff=p_cutoff,added_key=added_key)
        if pair_wise:
            # one verse one
            for a,b in itertools.combinations(df[group].map(str).unique(),2):
                df['Comparison']=df[group].apply(lambda x:"GroupA" if x==a else 'GroupB' if x==b else "Rest")
                cal_chi2_pvalue(adata,df.loc[df.Comparison.isin(['GroupA','GroupB'])],
                                cols=cols, n_group_a=n_count_dict[a],
                                n_group_b=n_count_dict[b],
                                comparison=a+" | "+b,
                                p_cutoff=p_cutoff,
                                added_key=added_key)

    adata.uns[added_key].sort_values('chi2_adjp',inplace=True)
    features=adata.uns[added_key].groupby('Comparison').apply(lambda x:x.head(topn).index.tolist()).sum()
    adata=adata[:,list(set(features))]
    return adata
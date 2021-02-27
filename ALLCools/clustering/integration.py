import numpy as np
import pandas as pd


def calculate_direct_confusion(left_part, right_part):
    """
    Given 2 dataframe for left/source and right/target dataset,
    calculate the direct confusion matrix based on co-cluster labels.
    Each dataframe only contain 2 columns, first is original cluster, second is co-cluster.
    The returned confusion matrix will be the form of source-cluster by target-cluster.

    Parameters
    ----------
    left_part
    right_part

    Returns
    -------

    """
    left_part = left_part.astype(str)
    original_left_name = left_part.columns[0]
    left_part.columns = ['cluster', 'co_cluster']

    right_part = right_part.astype(str)
    original_right_name = right_part.columns[0]
    if original_right_name == original_left_name:
        original_right_name += '_1'
    right_part.columns = ['cluster', 'co_cluster']

    left_confusion = left_part.groupby('cluster')['co_cluster'].value_counts().unstack()
    right_confusion = right_part.groupby('cluster')['co_cluster'].value_counts().unstack()

    left_confusion_portion = left_confusion.divide(left_confusion.sum(axis=1), axis=0).fillna(0)
    right_confusion_portion = right_confusion.divide(right_confusion.sum(axis=1), axis=0).fillna(0)

    records = []
    for left_cluster, left_row in left_confusion_portion.iterrows():
        for right_cluster, right_row in right_confusion_portion.iterrows():
            union_index = left_row.index | right_row.index
            left_row = left_row.reindex(union_index).fillna(0)
            right_row = right_row.reindex(union_index).fillna(0)
            overlap_value = np.min([left_row.values, right_row.values], axis=0).sum()
            records.append([left_cluster, right_cluster, overlap_value])

    flat_confusion_matrix = pd.DataFrame(records,
                                         columns=[original_left_name,
                                                  original_right_name,
                                                  'overlap_value'])
    confusion_matrix = flat_confusion_matrix.set_index([original_left_name, original_right_name]).squeeze().unstack()
    return confusion_matrix

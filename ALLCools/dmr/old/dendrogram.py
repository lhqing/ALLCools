import numpy as np
import pandas as pd


def extract_all_nodes(linkage, labels=None):
    """
    Given a linkage array output from scipy.cluster.hierarchy.linkage,
    calculate the left and right branches for all of the non-singleton nodes.
    """
    if isinstance(linkage, pd.DataFrame):
        linkage = linkage.values
    if isinstance(labels, pd.Series):
        labels = labels.tolist()

    cluster_dict = {}
    cur_cluster_id = len(linkage) + 1
    total_obs = len(linkage) + 1
    if labels is None:
        labels = list(range(total_obs))
    for left, right, distance, n_obs in linkage:
        left = int(left)
        right = int(right)
        n_obs = int(n_obs)

        cluster_dict[cur_cluster_id] = {'left': set(),
                                        'right': set()}
        if n_obs == 2:
            left = labels[left]
            right = labels[right]
            # merge of 2 original observations
            cluster_dict[cur_cluster_id]['left'].add(left)
            cluster_dict[cur_cluster_id]['right'].add(right)
        else:
            # left and/or right are cluster
            if left < total_obs:
                left = labels[left]
                cluster_dict[cur_cluster_id]['left'].add(left)
            else:
                # node are cluster
                cluster_dict[cur_cluster_id]['left'].update(cluster_dict[left]['left'])
                cluster_dict[cur_cluster_id]['left'].update(cluster_dict[left]['right'])
            if right < total_obs:
                right = labels[right]
                cluster_dict[cur_cluster_id]['right'].add(right)
            else:
                # node are cluster
                cluster_dict[cur_cluster_id]['right'].update(cluster_dict[right]['left'])
                cluster_dict[cur_cluster_id]['right'].update(cluster_dict[right]['right'])
        cur_cluster_id += 1
    return cluster_dict


def dmr_distance(a: np.ndarray, b: np.ndarray):
    overlap = (a == b)[(a * b) != 0]
    return 1 - (overlap.sum() / overlap.size)

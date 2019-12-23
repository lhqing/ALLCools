import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import pandas as pd
import numpy as np


def hclust_to_scipy_linkage(result, labels=None):
    """Turn R hclust result obj into scipy linkage matrix format"""
    # in hclust merge matrix, negative value is for singleton
    raw_linkage = pd.DataFrame(np.array(result[0]))
    nobs = raw_linkage.shape[0] + 1
    raw_linkage[2] = np.array(result[1])
    raw_linkage.index = raw_linkage.index + nobs

    # in hclust merge matrix, positive value is for non-singleton
    scipy_linkage = raw_linkage.copy()
    scipy_linkage[raw_linkage.iloc[:, :2] < 0] += nobs
    scipy_linkage[raw_linkage.iloc[:, :2] > 0] += (nobs - 1)
    total_obs = nobs

    # add the 4th col: number of singleton
    cluster_dict = {}
    if labels is None:
        labels = list(range(total_obs))
    for cur_cluster_id, (left, right, distance) in scipy_linkage.iterrows():
        left = int(left)
        right = int(right)
        cluster_dict[cur_cluster_id] = {'left': set(), 'right': set()}
        if (left < total_obs) and (right < total_obs):
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
                cluster_dict[cur_cluster_id]['left'].update(
                    cluster_dict[left]['left'])
                cluster_dict[cur_cluster_id]['left'].update(
                    cluster_dict[left]['right'])
            if right < total_obs:
                right = labels[right]
                cluster_dict[cur_cluster_id]['right'].add(right)
            else:
                # node are cluster
                cluster_dict[cur_cluster_id]['right'].update(
                    cluster_dict[right]['left'])
                cluster_dict[cur_cluster_id]['right'].update(
                    cluster_dict[right]['right'])
        cur_cluster_id += 1

    cluster_records = {}
    for cluster, _sub_dict in cluster_dict.items():
        total_n = len(_sub_dict['left']) + len(_sub_dict['right'])
        cluster_records[cluster] = total_n
    scipy_linkage[3] = pd.Series(cluster_records)
    return scipy_linkage


def pvclust(matrix_path, nboot=100, method_dist='correlation', method_hclust='average', cpu=10):
    base = importr("base")
    pvclust = importr("pvclust")
    dataframe = robjects.DataFrame.from_csvfile(matrix_path, sep=",", row_names=1)

    result = pvclust.pvclust(dataframe,
                             nboot=nboot,
                             method_dist=method_dist,
                             method_hclust=method_hclust,
                             parallel=cpu)

    edge_profile = result[1]
    edge_profile.to_csvfile(matrix_path + '.EdgeStats.csv')
    linkage_df = hclust_to_scipy_linkage(result[0])
    linkage_df.to_csv(matrix_path + '.Linkage.csv')
    return

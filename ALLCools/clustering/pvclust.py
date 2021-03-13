import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr, isinstalled
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
import joblib


def install_r_package(name):
    if not isinstalled(name):
        utils = importr('utils')
        utils.chooseCRANmirror(ind=1)
        utils.install_packages(StrVector([name]))


def _hclust_to_scipy_linkage(result, plot=True):
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

    # dendrogram
    orders = list(result[2])
    labels = list(result[3])
    # correct order of the final dendrogram
    r_order = [labels[i - 1] for i in orders]
    dendro = dendrogram(scipy_linkage.values, no_plot=True)
    python_order = pd.Series({a: b for a, b in zip(dendro['leaves'], r_order)}).sort_index().tolist()
    # python_order = [i[1:] for i in python_order]
    if plot:
        fig, ax = plt.subplots(dpi=300)
        dendro = dendrogram(scipy_linkage.values, labels=tuple(python_order), no_plot=False, ax=ax)
        ax.xaxis.set_tick_params(rotation=90)
    else:
        dendro = dendrogram(scipy_linkage.values, labels=tuple(python_order), no_plot=True)
    return scipy_linkage, python_order, dendro


class Dendrogram:
    def __init__(self,
                 nboot=1000,
                 method_dist='correlation',
                 method_hclust='average',
                 n_jobs=-1):
        self.nboot = nboot
        self.method_dist = method_dist
        self.method_hclust = method_hclust
        self.n_jobs = n_jobs

        self.linkage = None
        self.label_order = None
        self.dendrogram = None
        self.edge_stats = None

    def fit(self, data):
        """

        Parameters
        ----------
        data
            The data is in obs-by-var form, row is obs.

        Returns
        -------

        """
        importr("base")
        pvclust = importr("pvclust")
        with localconverter(ro.default_converter + pandas2ri.converter):
            # The data is in obs-by-var form, row is obs. Transpose for R.
            r_df = ro.conversion.py2rpy(data.T)
        if self.n_jobs == -1:
            self.n_jobs = True
        result = pvclust.pvclust(r_df,
                                 nboot=self.nboot,
                                 method_dist=self.method_dist,
                                 method_hclust=self.method_hclust,
                                 parallel=self.n_jobs)
        # dendrogram info
        hclust = result[0]
        linkage, label_order, dendro = _hclust_to_scipy_linkage(hclust, plot=False)
        self.linkage = linkage
        self.label_order = label_order
        self.dendrogram = dendro
        # scores of edges by pvclust bootstrap
        edge_stats = pd.DataFrame(result[1], index=result[1].colnames).T
        edge_stats.index = linkage.index
        child_dict = {}
        # pvclust edge stat is only for parents, here turn it into child basis
        for parent, (left, right, *_) in linkage.iterrows():
            child_dict[int(left)] = edge_stats.loc[parent]
            child_dict[int(right)] = edge_stats.loc[parent]
        self.edge_stats = pd.DataFrame(child_dict).T.sort_index()
        return

    def save(self, output_path):
        joblib.dump(self, output_path)

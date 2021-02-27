import pandas as pd
import numpy as np
from scipy.sparse import issparse
from pandas.api.types import is_categorical_dtype
from collections import defaultdict
import louvain
import scipy
import seaborn as sns
import patsy

from scanpy import settings as sett
from scanpy import logging as logg

import scanpy as sc
# for color
from scanpy.plotting.palettes import *

import matplotlib.pyplot as plt
from matplotlib import rcParams

# for reading/saving clf model
from sklearn.mixture import BayesianGaussianMixture
from sklearn import metrics
from sklearn.metrics.pairwise import euclidean_distances, pairwise_distances
from sklearn.model_selection import cross_val_score, cross_val_predict, cross_validate, train_test_split
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.svm import SVC

color_long = ['#e6194b', '#3cb44b', '#ffe119', '#0082c8', '#f58231', '#911eb4', \
              '#46f0f0', '#f032e6', '#d2f53c', '#fabebe', '#008080', '#e6beff', \
              '#aa6e28', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000080', \
              '#808080', '#000000', "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", \
              "#008941", "#006FA6", "#A30059", "#FFDBE5", "#7A4900", "#0000A6", \
              "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", \
              "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", \
              "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", \
              "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", \
              "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", \
              "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", \
              "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C", \
              "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81", \
              "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", \
              "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", \
              "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", \
              "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72"]


def run_BayesianGaussianMixture(Y, K):
    """
    For K-means clustering

    Input
    -----
    Y: the expression matrix
    K: number of clusters

    return
    -----
    clusters assigned to each cell.
    """
    gmm = BayesianGaussianMixture(K, max_iter=1000)
    gmm.fit(Y)
    return gmm.predict(Y)


def bhattacharyya_distance(repr1, repr2):
    """Calculates Bhattacharyya distance (https://en.wikipedia.org/wiki/Bhattacharyya_distance)."""
    sim = - np.log(np.sum([np.sqrt(p * q) for (p, q) in zip(repr1, repr2)]))
    assert not np.isnan(sim), 'Error: Similarity is nan.'
    if np.isinf(sim):
        # the similarity is -inf if no term in the review is in the vocabulary
        return 0
    return sim


def bhattacharyya_matrix(prob, flags=None):
    ndim = prob.shape[1]
    m = np.zeros((ndim, ndim))
    if flags is None:
        for i in np.arange(ndim):
            for j in np.arange(i + 1, ndim):
                val = -bhattacharyya_distance(prob[:, i], prob[:, j])
                m[i, j] = val
                m[j, i] = val
    else:
        for i in np.arange(ndim):
            for j in np.arange(i + 1, ndim):
                df = pd.DataFrame({'i': prob[:, i], 'j': prob[:, j], 'f': flags})
                df = df[df['f']]
                val = -bhattacharyya_distance(df['i'], df['j'])
                m[i, j] = val
                m[j, i] = val
    return m


def binary_accuracy(X, y, clf):
    y_pred = clf.predict(X)
    return y == y_pred


def normalize_confmat1(cmat, mod='1'):
    """
    Normalize the confusion matrix based on the total number of cells in each class
    x(i,j) = max(cmat(i,j)/diagnol(i),cmat(j,i)/diagnol(j))
    confusion rate between i and j is defined by the maximum ratio i is confused as j or j is confused as i.

    Input
    cmat: the confusion matrix

    return
    -----
    the normalized confusion matrix
    """
    dmat = cmat.values
    smat = np.diag(dmat)
    dim = cmat.shape[0]
    xmat = np.zeros([dim, dim])
    for i in range(dim):
        for j in range(i + 1, dim):
            if mod is '1':
                xmat[i, j] = xmat[j, i] = max(dmat[i, j] / smat[j], dmat[j, i] / smat[i])
            else:
                xmat[i, j] = xmat[j, i] = max(dmat[i, j] / smat[i], dmat[j, i] / smat[j])
    return xmat


def normalize_confmat2(cmat):
    """
    Normalize the confusion matrix based on the total number of cells.
    x(i,j) = max(cmat(i,j)+cmat(j,i)/N)
    N is total number of cells analyzed.
    Confusion rate between i and j is defined by the sum of i confused as j or j confused as i.
    Then divide by total number of cells.

    Input
    cmat: the confusion matrix

    return
    -----
    the normalized confusion matrix
    """
    emat = np.copy(cmat)
    s = np.sum(cmat)
    emat = emat + emat.T
    np.fill_diagonal(emat, 0)
    return emat * 1.0 / s


def cluster_adjmat(xmat,
                   resolution=1,
                   cutoff=0.1):
    """
    Cluster the groups based on the adjacent matrix.
    Use the cutoff to discretize the matrix used to construct the adjacent graph.
    Then cluster the graph using the louvain clustering with a resolution value.
    As the adjacent matrix is binary, the default resolution value is set to 1.

    Input
    -----
    xmat: `numpy.array` or sparse matrix
        the reference matrix/normalized confusion matrix
    cutoff: `float` optional (default: 0.1)
        threshold used to binarize the reference matrix
    resolution: `float` optional (default: 1.0)
        resolution parameter for louvain clustering

    return
    -----
    new group names.
    """
    g = sc._utils.get_igraph_from_adjacency((xmat > cutoff).astype(int), directed=False)
    print(g)
    part = louvain.find_partition(g, louvain.RBConfigurationVertexPartition,
                                  resolution_parameter=resolution)
    groups = np.array(part.membership)
    return groups


def msample(x, n, frac):
    """
    sample the matrix by number or by fraction.
    if the fraction is larger than the sample number, use number for sampling. Otherwise, use fraction.

    Input
    -----
    x: the matrix to be split
    n: number of vectors to be sampled
    frac: fraction of the total matrix to be sampled

    return
    -----
    sampled selection.
    """
    if len(x) <= np.floor(n / frac):
        if len(x) < 10: frac = 0.9
        return x.sample(frac=frac)
    else:
        return x.sample(n=n)


def train_test_split_per_type(X, y, n=100, frac=0.8):
    """
    This function is identical to train_test_split, but can split the data either based on number of cells or by fraction.

    Input
    -----
    X: `numpy.array` or sparse matrix
        the feature matrix
    y: `list of string/int`
        the class assignments
    n: `int` optional (default: 100)
        maximum number sampled in each label
    fraction: `float` optional (default: 0.8)
        Fraction of data included in the training set. 0.5 means use half of the data for training,
        if half of the data is fewer than maximum number of cells (n).

    return
    -----
    X_train, X_test, y_train, y_test
    """
    df = pd.DataFrame(y)
    df.index = np.arange(len(y))
    df.columns = ['class']
    c_idx = df.groupby('class').apply(lambda x: msample(x, n=n, frac=frac)).index.get_level_values(None)
    d_idx = ~np.isin(np.arange(len(y)), c_idx)
    return X[c_idx, :], X[d_idx, :], y[c_idx], y[d_idx]


# functions for SCCAF
def SCCAF_assessment(*args, **kwargs):
    """
    Assessment of clustering reliability using self-projection.
    It is the same as the self_projection function.
    """
    return self_projection(*args, **kwargs)


# need to check number of cells in each cluster of the training set.
def self_projection(X,
                    cell_types,
                    classifier="LR",
                    penalty='l1',
                    sparsity=0.5,
                    fraction=0.5,
                    solver='liblinear',
                    n=0,
                    cv=5,
                    whole=False,
                    n_jobs=None):
    # n = 100 should be good.
    """
    This is the core function for running self-projection.

    Input
    -----
    X: `numpy.array` or sparse matrix
        the expression matrix, e.g. ad.raw.X.
    cell_types: `list of String/int`
        the cell clustering assignment
    classifier: `String` optional (defatul: 'LR')
        a machine learning model in "LR" (logistic regression), \
        "RF" (Random Forest), "GNB"(Gaussion Naive Bayes), "SVM" (Support Vector Machine) and "DT"(Decision Tree).
    penalty: `String` optional (default: 'l2')
        the standardization mode of logistic regression. Use 'l1' or 'l2'.
    sparsity: `fload` optional (default: 0.5)
        The sparsity parameter (C in sklearn.linear_model.LogisticRegression) for the logistic regression model.
    fraction: `float` optional (default: 0.5)
        Fraction of data included in the training set. 0.5 means use half of the data for training,
        if half of the data is fewer than maximum number of cells (n).
    n: `int` optional (default: 100)
        Maximum number of cell included in the training set for each cluster of cells.
        only fraction is used to split the dataset if n is 0.
    cv: `int` optional (default: 5)
        fold for cross-validation on the training set.
        0 means no cross-validation.
    whole: `bool` optional (default: False)
        if measure the performance on the whole dataset (include training and test).
    n_jobs: `int` optional, number of threads to use with the different classifiers (default: None - unlimited).

    return
    -----
    y_prob, y_pred, y_test, clf
    y_prob: `matrix of float`
        prediction probability
    y_pred: `list of string/int`
        predicted clustering of the test set
    y_test: `list of string/int`
        real clustering of the test set
    clf: the classifier model.
    """
    # split the data into training and testing
    if n > 0:
        X_train, X_test, y_train, y_test = \
            train_test_split_per_type(X, cell_types, n=n, frac=(1 - fraction))
    else:
        X_train, X_test, y_train, y_test = \
            train_test_split(X, cell_types,
                             stratify=cell_types, test_size=fraction)  # fraction means test size
    # set the classifier
    if classifier == 'LR':
        clf = LogisticRegression(random_state=1, penalty=penalty, C=sparsity, multi_class="ovr", solver=solver)
    elif classifier == 'RF':
        clf = RandomForestClassifier(random_state=1, n_jobs=n_jobs)
    elif classifier == 'GNB':
        clf = GaussianNB()
    elif classifier == 'GPC':
        clf = GaussianProcessClassifier(n_jobs=n_jobs)
    elif classifier == 'SVM':
        clf = SVC(probability=True)
    elif classifier == 'SH':
        clf = SGDClassifier(loss='squared_hinge', n_jobs=n_jobs)
    elif classifier == 'PCP':
        clf = SGDClassifier(loss='perceptron', n_jobs=n_jobs)
    elif classifier == 'DT':
        clf = DecisionTreeClassifier()

    # mean cross validation score
    cvsm = 0
    if cv > 0:
        cvs = cross_val_score(clf, X_train, np.array(y_train), cv=cv, scoring='accuracy', n_jobs=n_jobs)
        cvsm = cvs.mean()
        print("Mean CV accuracy: %.4f" % cvsm)
    # accuracy on cross validation and on test set
    clf.fit(X_train, y_train)
    accuracy = clf.score(X_train, y_train)
    print("Accuracy on the training set: %.4f" % accuracy)
    accuracy_test = clf.score(X_test, y_test)
    print("Accuracy on the hold-out set: %.4f" % accuracy_test)

    # accuracy of the whole dataset
    if whole:
        accuracy = clf.score(X, cell_types)
        print("Accuracy on the whole set: %.4f" % accuracy)

    # get predicted probability on the test set
    y_prob = None
    if not classifier in ['SH', 'PCP']:
        y_prob = clf.predict_proba(X_test)
    y_pred = clf.predict(X_test)

    return y_prob, y_pred, y_test, clf, cvsm, accuracy_test


def make_unique(dup_list):
    """
    Make a name list unique by adding suffix "_%d". This function is identical to the make.unique function in R.

    Input
    -----
    dup_list: a list

    return
    -----
    a unique list with the same length as the input.
    """
    from collections import Counter

    counter = Counter()
    deduped = []
    for name in dup_list:
        new = str(name) + "_%s" % str(counter[name]) if counter[name] else name
        counter.update({name: 1})
        deduped.append(new)
    return deduped


def confusion_matrix(y_test, y_pred, clf, labels=None):
    """
    Get confusion matrix based on the test set.

    Input
    -----
    y_test, y_pred, clf: same as in self_projection

    return
    -----
    the confusion matrix
    """
    if labels is None: labels = clf.classes_
    df = pd.DataFrame.from_records(metrics.confusion_matrix(y_test, y_pred, labels=labels))
    df.index = labels
    df.columns = labels
    df.index.name = 'Real'
    df.columns.name = 'Predicted'
    return df


def per_cluster_accuracy(mtx, ad=None, clstr_name='louvain'):
    """
    Measure the accuracy of each cluster and put into a metadata slot.
    So the reliability of each cluster can be visualized.

    Input
    -----
    mtx: `pandas.dataframe`
        the confusion matrix
    ad: `AnnData`
        anndata object
    clstr_name: `String`
        the name of the clustering

    return
    -----
    """
    ds = None
    if not ad is None:
        ds = (np.diag(mtx.values) / mtx.sum(0)).fillna(0)
        rel_name = '%s_reliability' % clstr_name
        ad.obs[rel_name] = ad.obs[clstr_name].astype('category')
        ad.obs[rel_name].cat.categories = make_unique(ds)
        ad.obs[rel_name] = ad.obs[rel_name].astype(str).str.split("_").str[0]
        ad.obs[rel_name] = ad.obs[rel_name].astype(float)
    return ds


def per_cell_accuracy(X, cell_types, clf):
    y_prob = clf.predict_proba(X)
    df = pd.DataFrame(y_prob, index=cell_types, columns=clf.classes_).sort_index().T
    df.index.name = 'Predicted'
    dy = np.empty([0])
    for cell in df.index:
        x = np.array(df.loc[cell][cell].values)
        dy = np.concatenate((dy, x))
    return dy / np.max(df.values, 0)


def get_topmarkers(clf, names, topn=10):
    """
    Get the top weighted features from the logistic regressioin model.

    Input
    -----
    clf: the logistic regression classifier
    names: `list of Strings`
        the names of the features (the gene names).
    topn: `int`
        number of top weighted featured to be returned.

    return
    -----
    list of markers for each of the cluster.
    """
    marker_genes = pd.DataFrame({
        'cell_type': clf.classes_[clf.coef_.argmax(0)],
        'gene': names,
        'weight': clf.coef_.max(0)
    })

    top_markers = \
        marker_genes \
            .query('weight > 0.') \
            .sort_values('weight', ascending=False) \
            .groupby('cell_type') \
            .head(topn) \
            .sort_values(['cell_type', 'weight'], ascending=[True, False])
    return top_markers


def eu_distance(X, gp1, gp2, cell):
    """
    Measure the euclidean distance between two groups of cells and the third group.

    Input
    -----
    X: `np.array` or `sparse matrix`
        the total expression matrix
    gp1: `bool list`
        group1 of cells
    gp2: `bool list`
        group2 of cells
    cell: `bool list`
        group3 of cells, the group to be compared with gp1 and gp2.

    return
    -----
    `float value`
    the average distance difference.
    """
    d1 = euclidean_distances(X[gp1 & (~cell), :], X[cell, :])
    d2 = euclidean_distances(X[gp2 & (~cell), :], X[cell, :])
    df1 = pd.DataFrame(d1[:, 0], columns=['distance'])
    df1['type'] = 'cell'
    df2 = pd.DataFrame(d2[:, 0], columns=['distance'])
    df2['type'] = 'cell_pred'
    df = pd.concat([df1, df2])
    m1 = d1.mean()
    m2 = d2.mean()
    print('%f - %f' % (m1, m2))
    return df


def get_distance_matrix(X, clusters, labels=None, metric='euclidean'):
    """
    Get the mean distance matrix between all clusters.

    Input
    -----
    X: `np.array` or `sparse matrix`
        the total expression matrix
    clusters: `string list`
        the assignment of the clusters
    labels: `string list`
        the unique labels of the clusters
    metric: `string` (optional, default: euclidean)
        distance metrics, see (http://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise_distances.html)

    return
    -----
    `np.array`
        the all-cluster to all-cluster distance matrix.
    """
    if labels is None:
        labels = np.unique(clusters)
    centers = []
    if scipy.sparse.issparse(X):
        for cl in labels:
            centers.append(np.array(X[np.where(clusters == cl)[0], :].mean(0))[0, :])
    else:
        for cl in labels:
            centers.append(np.array(X[np.where(clusters == cl)[0], :].mean(0)))
    return pairwise_distances(np.array(centers), metric=metric)


def merge_cluster(ad, old_id, new_id, groups):
    ad.obs[new_id] = ad.obs[old_id]
    ad.obs[new_id] = ad.obs[new_id].astype('category')
    ad.obs[new_id].cat.categories = make_unique(groups.astype(str))
    ad.obs[new_id] = ad.obs[new_id].str.split('_').str[0]
    return ad


def find_high_resolution(ad, resolution=4, n=100):
    cut = resolution
    while cut > 0.5:
        print("clustering with resolution: %.1f" % cut)
        sc.tl.leiden(ad, resolution=cut)
        ad.obs['leiden_res%.1f' % cut] = ad.obs['leiden']
        if ad.obs['leiden'].value_counts().min() > n:
            break
        cut -= 0.5


def get_connection_matrix(ad_obs, key1, key2):
    df = ad_obs.groupby([key1, key2]).size().to_frame().reset_index()
    df.columns = [key1, key2, 'counts']
    df2 = ad_obs[key2].value_counts()
    df['size'] = df2[df[key2].tolist()].tolist()
    df['percent'] = df['counts'] / df['size']
    df = df[df['percent'] > 0.1]
    df2 = df.groupby(key1).size()
    df2 = df2[df2 > 1]
    df = df[df[key1].isin(df2.index)]

    dim = len(ad_obs[key2].unique())
    mat = pd.DataFrame(np.zeros([dim, dim])).astype(int)
    mat.columns = mat.columns.astype(str)
    mat.index = mat.index.astype(str)

    import itertools
    grouped = df.groupby(key1)
    for name, group in grouped:
        x = group[key2]
        if len(x) > 0:
            for i, j in list(itertools.combinations(x.tolist(), 2)):
                mat.loc[i, j] = mat.loc[j, i] = 1
    return mat


def SCCAF_optimize_all(adata,
                       start_groups=None,
                       min_acc=0.9,
                       r1_norm_cutoff=0.5,
                       r2_norm_cutoff=0.05,
                       R1norm_step=0.01,
                       R2norm_step=0.001,
                       min_iter=3,
                       max_iter=10,
                       *args, **kwargs):
    """
    adata: `AnnData`
        The AnnData object of the expression profile.
    min_acc: `float` optional (default: 0.9)
        The minimum self-projection accuracy to be optimized for.
        e.g., 0.9 means the clustering optimization (merging process)
        will not stop until the self-projection accuracy is above 90%.
    R1norm_cutoff: `float` optional (default: 0.5)
        The start cutoff for R1norm of confusion matrix.
        e.g., 0.5 means use 0.5 as a cutoff to discretize the confusion matrix after R1norm.
        the discretized matrix is used to construct the connection graph for clustering optimization.
    R2norm_cutoff: `float` optional (default: 0.05)
        The start cutoff for R2norm of confusion matrix.
    R1norm_step: `float` optional (default: 0.01)
        The reduce step for minimum R1norm value.
        Each round of optimization calls the function `SCCAF_optimize`.
        The start of the next round of optimization is based on a new
        cutoff for the R1norm of the confusion matrix. This cutoff is
        determined by the minimum R1norm value in the previous round minus the R1norm_step value.
    R2norm_step: `float` optional (default: 0.001)
        The reduce step for minimum R2norm value.
    """
    prefix = start_groups
    cur_iter = 0
    if not (start_groups in adata.obs.keys()):
        raise ValueError(f"`adata.obs['{start_groups}']` doesn't exist. Please assign the initial clustering first.")
    adata.obs[f'{prefix}_Round{cur_iter}'] = adata.obs[start_groups]

    old_n_cluster = len(adata.obs[f'{prefix}_Round{cur_iter}'].unique())
    # 'while acc < min_acc:
    for i in range(1, max_iter + 1):
        print(f"Current iteration: {i}")
        print(f"R1norm_cutoff: {r1_norm_cutoff:.3f}")
        print(f"R2norm_cutoff: {r2_norm_cutoff:.3f}")
        print("======================")
        adata, m1, m2, acc, start_iter = SCCAF_optimize(ad=adata,
                                                        r1_norm_cutoff=r1_norm_cutoff,
                                                        r2_norm_cutoff=r2_norm_cutoff,
                                                        start_iter=start_iter,
                                                        min_acc=min_acc,
                                                        prefix=prefix,
                                                        *args, **kwargs)
        print(f"m1: {m1:.3f}")
        print(f"m2: {m2:.3f}")
        print(f"Accuracy: {acc:.3f}")
        r1_norm_cutoff = m1 - R1norm_step
        r2_norm_cutoff = m2 - R2norm_step

        new_n_cluster = len(adata.obs['%s_result' % prefix].unique())
        if new_n_cluster >= old_n_cluster and i >= min_iter:
            print("Cluster merging converged.")
            break

        if acc >= min_acc:
            print(f'Accuracy passed the threshold {min_acc:.3f}.')
            break
    return adata


def SCCAF_optimize(ad,
                   prefix='L1',
                   use='raw',
                   use_projection=False,
                   R1norm_only=False,
                   R2norm_only=False,
                   dist_only=False,
                   dist_not=True,
                   plot=True,
                   basis='umap',
                   plot_dist=False,
                   plot_cmat=False,
                   mod='1',
                   low_res=None,
                   c_iter=3,
                   n_iter=10,
                   n_jobs=None,
                   start_iter=0,
                   sparsity=0.5,
                   n=100,
                   fraction=0.5,
                   r1_norm_cutoff=0.1,
                   r2_norm_cutoff=1,
                   dist_cutoff=8,
                   classifier="LR",
                   mplotlib_backend=None,
                   min_acc=1):
    """
    This is a self-projection confusion matrix directed cluster optimization function.

    Input
    -----
    ad: `AnnData`
        The AnnData object of the expression profile.
    prefix: `String`, optional (default: 'L1')
        The name of the optimization, which set as a prefix.
        e.g., the prefix = 'L1', the start round of optimization clustering is based on
        'L1_Round0'. So we need to assign an over-clustering state as a start point.
        e.g., ad.obs['L1_Round0'] = ad.obs['louvain']
    use: `String`, optional (default: 'raw')
        Use what features to train the classifier. Three choices:
        'raw' uses all the features;
        'hvg' uses the highly variable genes in the anndata object ad.var_names slot;
        'pca' uses the PCA data in the anndata object ad.obsm['X_pca'] slot.
    R1norm_only: `bool` optional (default: False)
        If only use the confusion matrix(R1norm) for clustering optimization.
    R2norm_only: `bool` optional (default: False)
        If only use the confusion matrix(R2norm) for clustering optimization.
    dist_only: `bool` optional (default: False)
        If only use the distance matrix for clustering optimization.
    dist_not: `bool` optional (default: True)
        If not use the distance matrix for clustering optimization.
    plot: `bool` optional (default: True)
        If plot the self-projectioin results, ROC curves and confusion matrices,
        during the optimization.
    plot_tsne: `bool` optional (default: False)
        If plot the self-projectioin results as tSNE. If False, the results are plotted as UMAP.
    plot_dist: `bool` optional (default: False)
        If make a scatter plot of the distance compared with the confusion rate for each of the cluster.
    plot_cmat: `bool` optional (default: False)
        plot the confusion matrix or not.
    mod: `string` optional (default: '1')
        two directions of normalization of confusion matrix for R1norm.
    c_iter: `int` optional (default: 3)
        Number of iterations of sampling for the confusion matrix.
        The minimum value of confusion rate in all the iterations is used as the confusion rate between two clusters.
    n_iter： `int` optional (default: 10)
        Maximum number of iterations(Rounds) for the clustering optimization.
    start_iter： `int` optional (default: 0)
        The start round of the optimization. e.g., start_iter = 3,
        the optimization will start from ad.obs['%s_3'%prefix].
    sparsity: `fload` optional (default: 0.5)
        The sparsity parameter (C in sklearn.linear_model.LogisticRegression) for the logistic regression model.
    n: `int` optional (default: 100)
        Maximum number of cell included in the training set for each cluster of cells.
    n_jobs: `int` number of jobs/threads to use (default: None - unlimited).
    fraction: `float` optional (default: 0.5)
        Fraction of data included in the training set. 0.5 means use half of the data for training,
        if half of the data is fewer than maximum number of cells (n).
    R1norm_cutoff: `float` optional (default: 0.1)
        The cutoff for the confusion rate (R1norm) between two clusters.
        0.1 means we allow maximum 10% of the one cluster confused as another cluster.
    R2norm_cutoff: `float` optional (default: 1.0)
        The cutoff for the confusion rate (R2norm) between two clusters.
        1.0 means the confusion between any two cluster should not exceed 1% of the total number of cells.
    dist_cutoff: `float` optional (default: 8.0)
        The cutoff for the euclidean distance between two clusters of cells.
        8.0 means the euclidean distance between two cell types should be greater than 8.0.
    low_res: `str` optional
		the clustering boundary for under-clustering. Set a low resolution in louvain/leiden clustering and give
		the key as the underclustering boundary.
    classifier: `String` optional (default: 'LR')
        a machine learning model in "LR" (logistic regression), \
        "RF" (Random Forest), "GNB"(Gaussion Naive Bayes), "SVM" (Support Vector Machine) and "DT"(Decision Tree).
    mplotlib_backend: `matplotlib.backends.backend_pdf` optional
        MatPlotLib multi-page backend object instance, previously initialised (currently the only type supported is
        PdfPages).
    min_acc: `float`
		the minimum total accuracy to be achieved. Above this threshold, the optimization will stop.

    return
    -----
    The modified anndata object, with a slot "%s_result"%prefix
        assigned as the clustering optimization results.
    """
    # the space to use
    X = None
    if use == 'raw':
        X = ad.raw.X
    elif use == 'pca':
        if 'X_pca' not in ad.obsm.keys():
            raise ValueError("`adata.obsm['X_pca']` doesn't exist. Run `sc.pp.pca` first.")
        X = ad.obsm['X_pca']
    else:
        X = ad[:, ad.var['highly_variable']].X

    for i in range(start_iter, start_iter + n_iter):
        print("Round%d ..." % (i + 1))
        old_id = '%s_Round%d' % (prefix, i)
        new_id = '%s_Round%d' % (prefix, i + 1)

        labels = np.sort(ad.obs[old_id].unique().astype(int)).astype(str)

        # optimize
        y_prob, y_pred, y_test, clf, cvsm, acc = self_projection(
            X, ad.obs[old_id], sparsity=sparsity, n=n,
            fraction=fraction, classifier=classifier, n_jobs=n_jobs)
        accs = [acc]
        ad.obs['%s_self-projection' % old_id] = clf.predict(X)

        cmat = confusion_matrix(y_test, y_pred, clf, labels=labels)
        xmat = normalize_confmat1(cmat, mod)
        xmats = [xmat]
        cmats = [np.array(cmat)]
        old_id1 = old_id
        if use_projection:
            old_id1 = '%s_self-projection' % old_id
        for j in range(c_iter - 1):
            y_prob, y_pred, y_test, clf, _, acc = self_projection(
                X, ad.obs[old_id1], sparsity=sparsity, n=n,
                fraction=fraction, classifier=classifier, cv=0,
                n_jobs=n_jobs)
            accs.append(acc)
            cmat = confusion_matrix(y_test, y_pred, clf, labels=labels)
            xmat = normalize_confmat1(cmat, mod)
            xmats.append(xmat)
            cmats.append(np.array(cmat))
        R1mat = np.minimum.reduce(xmats)
        R2mat = normalize_confmat2(np.minimum.reduce(cmats))

        m1 = np.max(R1mat)
        if np.isnan(m1):
            m1 = 1.
        m2 = np.max(R2mat)
        print("Max R1mat: %f" % m1)
        print("Max R2mat: %f" % m2)

        if np.min(accs) > min_acc:
            ad.obs['%s_result' % prefix] = ad.obs[old_id]
            print("Converge SCCAF_optimize min_acc!")
            break
        print("min_acc: %f" % np.min(accs))

        if R1norm_only:
            groups = cluster_adjmat(R1mat, cutoff=r1_norm_cutoff)
        elif R2norm_only:
            groups = cluster_adjmat(R2mat, cutoff=r2_norm_cutoff)
        else:
            if not low_res is None:
                conn_mat = get_connection_matrix(ad_obs=ad.obs, key1=low_res, key2=old_id)
                zmat = np.minimum.reduce([(R1mat > r1_norm_cutoff), conn_mat.values])
                groups = cluster_adjmat(zmat, cutoff=0)
            else:
                zmat = np.maximum.reduce([(R1mat > r1_norm_cutoff), (R2mat > r2_norm_cutoff)])
                groups = cluster_adjmat(zmat, cutoff=0)

        if len(np.unique(groups)) == len(ad.obs[old_id].unique()):
            ad.obs['%s_result' % prefix] = ad.obs[old_id]
            print("Converge SCCAF_optimize no. cluster!")
            break

        merge_cluster(ad, old_id1, new_id, groups)

        if len(np.unique(groups)) <= 1:
            ad.obs['%s_result' % prefix] = ad.obs[new_id]
            print("no clustering!")
            break

    return ad, m1, m2, np.min(accs), i

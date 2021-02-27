from concurrent.futures import ProcessPoolExecutor, as_completed

import anndata
import joblib
import leidenalg
import numpy as np
import pandas as pd
from hdbscan import HDBSCAN
from imblearn.ensemble import BalancedRandomForestClassifier
from natsort import natsorted
from scanpy.neighbors import Neighbors
from sklearn.metrics import confusion_matrix, pairwise_distances, balanced_accuracy_score
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import cross_val_predict


def _r1_normalize(cmat):
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
    dmat = cmat
    smat = np.diag(dmat) + 1  # in case some label has no correct prediction (0 in diag)
    dim = cmat.shape[0]
    xmat = np.zeros([dim, dim])
    for i in range(dim):
        for j in range(i + 1, dim):
            xmat[i, j] = xmat[j, i] = max(dmat[i, j] / smat[j], dmat[j, i] / smat[i])
    return xmat


def _r2_normalize(cmat):
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


def _leiden_runner(g, random_states, partition_type, **partition_kwargs):
    """run leiden clustering len(random_states) times with different random states,
    return all clusters as a pd.DataFrame"""
    results = []
    for seed in random_states:
        part = leidenalg.find_partition(g, partition_type, seed=seed, **partition_kwargs)
        groups = np.array(part.membership)
        groups = pd.Categorical(
            values=groups.astype('U'),
            categories=natsorted(np.unique(groups).astype('U')))
        results.append(groups)
    result_df = pd.DataFrame(results, columns=random_states)
    return result_df


def _split_train_test_per_group(x, y, frac, max_train, random_state):
    y_series = pd.Series(y)
    # split train test per group
    train_idx = []
    test_idx = []
    outlier_idx = []
    for cluster, sub_series in y_series.groupby(y_series):
        if (cluster == -1) or (sub_series.size < 3):
            outlier_idx += sub_series.index.tolist()
        else:
            n_train = max(1, min(max_train, int(sub_series.size * frac)))
            is_train = sub_series.index.isin(sub_series.sample(n_train, random_state=random_state).index)
            train_idx += sub_series.index[is_train].tolist()
            test_idx += sub_series.index[~is_train].tolist()
    x_train = x[train_idx]
    y_train = y[train_idx]
    x_test = x[test_idx]
    y_test = y[test_idx]
    return x_train, y_train, x_test, y_test


def _single_supervise_evaluation(
        clf,
        x_train,
        y_train,
        x_test,
        y_test,
        r1_norm_step=0.01,
        r2_norm_step=0.001,
        plot=True):
    # fit model
    clf.fit(x_train, y_train)

    # calc accuracy
    y_train_pred = clf.predict(x_train)
    accuracy_train = balanced_accuracy_score(y_true=y_train, y_pred=y_train_pred)
    print(f"Balanced accuracy on the training set: {accuracy_train:.3f}")
    y_test_pred = clf.predict(x_test)
    accuracy_test = balanced_accuracy_score(y_true=y_test, y_pred=y_test_pred)
    print(f"Balanced accuracy on the hold-out set: {accuracy_test:.3f}")

    # get confusion matrix
    y_pred = clf.predict(x_test)
    cmat = confusion_matrix(y_test, y_pred)

    # normalize confusion matrix
    r1_cmat = _r1_normalize(cmat)
    r2_cmat = _r2_normalize(cmat)
    m1 = np.max(r1_cmat)
    if np.isnan(m1):
        m1 = 1.
    m2 = np.max(r2_cmat)
    r1_norm_cutoff = m1 - r1_norm_step
    r2_norm_cutoff = m2 - r2_norm_step

    # final binary matrix to calculate which clusters need to be merged
    cluster_map = {}
    judge = np.maximum.reduce([(r1_cmat > r1_norm_cutoff),
                               (r2_cmat > r2_norm_cutoff)])
    if judge.sum() > 0:
        rows, cols = np.where(judge)
        edges = zip(rows.tolist(), cols.tolist())
        g = nx.Graph()
        g.add_edges_from(edges)
        for comp in nx.connected_components(g):
            to_label = comp.pop()
            for remain in comp:
                cluster_map[remain] = to_label

    if plot:
        fig, axes = plt.subplots(figsize=(12, 3), dpi=300, ncols=4)
        ax = axes[0]
        sns.heatmap(np.log1p(cmat), ax=ax)
        ax.set_title(f'confusion matrix')
        ax = axes[1]
        sns.heatmap(r1_cmat, ax=ax)
        ax.set_title(f'R1 Norm. cutoff {r1_norm_cutoff:.3f}')
        ax = axes[2]
        sns.heatmap(r2_cmat, ax=ax)
        ax.set_title(f'R2 Norm. cutoff {r2_norm_cutoff:.3f}')
        ax = axes[3]
        sns.heatmap(judge, ax=ax)
        ax.set_title(f'Merge')
        fig.tight_layout()
        fig.show()
    return clf, accuracy_test, cluster_map


class ConsensusClustering:
    def __init__(self,
                 model=None,
                 n_neighbors=25,
                 metric='euclidean',
                 min_cluster_size=10,
                 consensus_rate=0.9,
                 repeats=200,
                 resolution=1,
                 random_state=0,
                 train_frac=0.5,
                 train_max_n=500,
                 max_iter=10,
                 target_score=0.95,
                 n_jobs=-1,
                 outlier_cell_cutoff=0.25,
                 outlier_label='Outlier',
                 plot=True):
        # input metrics
        self.min_cluster_size = min_cluster_size
        self.consensus_rate = consensus_rate
        self.leiden_repeats = repeats
        self.leiden_resolution = resolution
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.n_neighbors = n_neighbors
        self.knn_metric = metric
        self.train_frac = train_frac
        self.train_max_n = train_max_n
        self.plot = plot
        self.max_iter = max_iter
        self.target_score = target_score
        self.outlier_cell_cutoff = outlier_cell_cutoff
        self.outlier_label = outlier_label

        self.n_obs, self.n_pcs = None, None
        self.X = None
        self._neighbors = None

        # multiple leiden clustering
        self.leiden_result_df = None
        self._multi_leiden_clusters = None

        # model training and outlier rescue
        self.supervise_model = model
        self._label_with_leiden_outliers = None
        self.label = None
        self.label_proba = None
        self.final_accuracy = None
        return

    def fit_predict(self, x, leiden_kwds=None):
        self.n_obs, self.n_pcs = x.shape
        self.X = x

        # Construct KNN graph
        print('Computing nearest neighbor graph')
        self.compute_neighbors()

        # repeat Leiden clustering with different random seeds
        print('Computing multiple clustering with different random seeds')
        kwds = {}
        if leiden_kwds is None:
            leiden_kwds = kwds
        else:
            leiden_kwds.update(kwds)
        self.multi_leiden_clustering(**leiden_kwds)

        # merge the over clustering version by supervised learning
        self.supervise_learning()

        # assign outliers
        self.final_evaluation()

    def compute_neighbors(self):
        # nearest neighbors graph
        adata = anndata.AnnData(X=None,
                                obs=pd.DataFrame([], index=[f'obs{i}' for i in range(self.n_obs)]),
                                var=pd.DataFrame([], index=[f'var{i}' for i in range(self.n_pcs)]))
        adata.obsm['X_pca'] = self.X
        # here neighbors should only use PCs
        self._neighbors = Neighbors(adata=adata)
        self._neighbors.compute_neighbors(
            n_neighbors=self.n_neighbors,
            knn=True,
            n_pcs=self.n_pcs,
            use_rep='X_pca',
            method='umap',
            metric=self.knn_metric,
            random_state=self.random_state)
        return

    def multi_leiden_clustering(self,
                                partition_type=None,
                                partition_kwargs=None,
                                use_weights=True,
                                n_iterations=-1):
        """Modified from scanpy"""
        if self._neighbors is None:
            raise ValueError('Run compute_neighbors first before multi_leiden_clustering')

        # convert neighbors to igraph
        g = self._neighbors.to_igraph()

        # generate n different seeds for each single leiden partition
        np.random.seed(self.random_state)
        leiden_repeats = self.leiden_repeats
        n_jobs = self.n_jobs
        random_states = np.random.choice(range(99999999), size=leiden_repeats, replace=False)
        step = max(int(leiden_repeats / n_jobs), 20)
        random_state_chunks = [random_states[i: min(i + step, leiden_repeats)]
                               for i in range(0, leiden_repeats, step)]

        results = []
        print(f'Repeating leiden clustering {leiden_repeats} times')
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            future_dict = {}
            for i, random_state_chunk in enumerate(random_state_chunks):
                # flip to the default partition type if not over writen by the user
                if partition_type is None:
                    partition_type = leidenalg.RBConfigurationVertexPartition
                # prepare find_partition arguments as a dictionary, appending to whatever the user provided
                # it needs to be this way as this allows for the accounting of a None resolution
                # (in the case of a partition variant that doesn't take it on input)
                if partition_kwargs is None:
                    partition_kwargs = {}
                else:
                    if 'seed' in partition_kwargs:
                        print('Warning: seed in the partition_kwargs will be ignored, use seed instead.')
                        del partition_kwargs['seed']
                if use_weights:
                    partition_kwargs['weights'] = np.array(g.es['weight']).astype(np.float64)
                partition_kwargs['n_iterations'] = n_iterations
                partition_kwargs['resolution_parameter'] = self.leiden_resolution
                # clustering proper
                future = executor.submit(_leiden_runner,
                                         g=g,
                                         random_states=random_state_chunk,
                                         partition_type=partition_type,
                                         **partition_kwargs)
                future_dict[future] = random_state_chunks

            for future in as_completed(future_dict):
                _ = future_dict[future]
                try:
                    data = future.result()
                    results.append(data)
                except Exception as exc:
                    print(f'_leiden_runner generated an exception: {exc}')
                    raise exc
        total_result = pd.concat(results, axis=1, sort=True)
        self.leiden_result_df = total_result
        cluster_count = self.leiden_result_df.apply(lambda i: i.unique().size)
        print(f'Found {cluster_count.min()} - {cluster_count.max()} clusters, '
              f'mean {cluster_count.mean():.1f}, std {cluster_count.std():.2f}')

        # create a over-clustering version based on all the leiden runs
        print('Summarizing multiple clustering results')
        self.summarize_multi_leiden()
        return

    def summarize_multi_leiden(self):
        # here we don't rely on hdbscan clustering, just use it to reduce the pairwise distances calculation
        # the resulting clusters are simply disconnected components based on hamming dist graph
        # This results an over-clustering, which will be merged by the supervise learning step.
        hdbscan = HDBSCAN(min_cluster_size=2,
                          min_samples=1,
                          metric='hamming',
                          cluster_selection_method='eom',
                          allow_single_cluster=True)
        hdbscan.fit(self.leiden_result_df)
        clusters = {}
        cur_num = 0
        for _, sub_df in self.leiden_result_df.groupby(
                pd.Series(hdbscan.labels_, index=self.leiden_result_df.index)
        ):
            pairwise_dist = pairwise_distances(sub_df, metric='hamming')
            # create a graph, cells within hamming_dist_cutoff are connected
            rows, cols = np.where(pairwise_dist < self.consensus_rate)
            edges = zip(sub_df.index[rows].tolist(), sub_df.index[cols].tolist())
            g = nx.Graph()
            g.add_edges_from(edges)
            for comp in nx.connected_components(g):
                if len(comp) >= self.min_cluster_size:
                    for node in comp:
                        # each disconnected component assigned to a cluster
                        clusters[node] = cur_num
                    cur_num += 1
                else:
                    for node in comp:
                        clusters[node] = -1

        clusters = pd.Series(clusters).sort_index()
        print(f'{(clusters != -1).sum()} cells assigned to {clusters.unique().size} raw clusters')
        print(f'{(clusters == -1).sum()} cells are outliers')
        self._multi_leiden_clusters = clusters.values
        return

    def _create_model(self,
                      n_estimators=1000):
        clf = BalancedRandomForestClassifier(
            n_estimators=n_estimators,
            criterion='gini',
            max_depth=None,
            min_samples_split=2,
            min_samples_leaf=2,
            min_weight_fraction_leaf=0.0,
            max_features='auto',
            max_leaf_nodes=None,
            min_impurity_decrease=0.0,
            bootstrap=True,
            oob_score=False,
            sampling_strategy='auto',
            replacement=False,
            n_jobs=self.n_jobs,
            random_state=self.random_state,
            verbose=0,
            warm_start=False,
            class_weight=None)
        return clf

    def supervise_learning(self):
        if self._multi_leiden_clusters is None:
            raise ValueError('Run multi_leiden_clustering first to get a '
                             'clustering assignment before run supervise_learning.')

        n_cluster = np.unique(self._multi_leiden_clusters[self._multi_leiden_clusters != -1]).size
        if n_cluster == 1:
            print('There is only one cluster except for outliers, can not train supervise model on that.')
            return
        print(f'\n=== Start supervise model training and cluster merging ===')

        x = self.X
        cur_y = self._multi_leiden_clusters.copy()
        score = None

        if self.supervise_model is None:
            # create default model if no model provided
            clf = self._create_model(n_estimators=500)
        else:
            clf = self.supervise_model
        for cur_iter in range(1, self.max_iter + 1):
            print(f'\n=== iteration {cur_iter} ===')
            n_labels = np.unique(cur_y[cur_y != -1]).size
            print(f'{n_labels} non-outlier labels')
            if n_labels < 2:
                print(f'Stop iteration because only {n_labels} cluster remain.')

            x_train, y_train, x_test, y_test = _split_train_test_per_group(
                x=x,
                y=cur_y,
                frac=self.train_frac,
                max_train=self.train_max_n,
                random_state=self.random_state + cur_iter  # every time train-test split got a different random state
            )
            clf, score, cluster_map = _single_supervise_evaluation(clf,
                                                                   x_train,
                                                                   y_train,
                                                                   x_test,
                                                                   y_test,
                                                                   r1_norm_step=0.01,
                                                                   r2_norm_step=0.001)
            if self.target_score < score:
                print('Stop iteration because score > target score.')
                break
            if len(cluster_map) > 0:
                print(f'Merging {len(cluster_map)} clusters.')
                cur_y = pd.Series(cur_y).apply(lambda i: cluster_map[i] if i in cluster_map else i)
                # renumber labels from large to small
                ordered_map = {c: i for i, c in enumerate(cur_y[cur_y != -1].value_counts().index)}
                cur_y = pd.Series(cur_y).apply(lambda i: ordered_map[i] if i in ordered_map else i).values
            else:
                print('Stop iteration because there is no cluster to merge')
                break
        else:
            print('Stop iteration because reaching maximum iteration.')
        self._label_with_leiden_outliers = cur_y
        self.label = cur_y
        self.supervise_model = clf
        self.final_accuracy = score
        return

    def final_evaluation(self):
        print(f'\n=== Assign final labels ===')
        # final evaluation of non-outliers using cross val predict
        _x = self.X[self.label != -1]
        _idx = np.where(self.label != -1)[0]
        _y = self.label[self.label != -1]
        final_predict_proba = cross_val_predict(self.supervise_model,
                                                _x,
                                                y=_y,
                                                method='predict_proba',
                                                n_jobs=self.n_jobs,
                                                verbose=0,
                                                cv=10)
        final_predict = pd.Series(np.argmax(final_predict_proba, axis=1), index=_idx)
        final_cell_proba = pd.Series(np.max(final_predict_proba, axis=1), index=_idx)
        final_acc = balanced_accuracy_score(_y, final_predict.values)
        print(f'Final CV Accuracy on non-outliers: {final_acc:.3f}')

        # predict outliers
        outlier_x = self.X[self.label == -1]
        outlier_idx = np.where(self.label == -1)[0]
        if len(outlier_idx) != 0:
            outlier_predict = pd.Series(self.supervise_model.predict(outlier_x),
                                        index=outlier_idx)
            outlier_proba = self.supervise_model.predict_proba(outlier_x)
            outlier_cell_proba = pd.Series(outlier_proba.max(axis=1),
                                           index=outlier_idx)
            final_predict = pd.concat([final_predict, outlier_predict]).sort_index()
            final_cell_proba = pd.concat([final_cell_proba, outlier_cell_proba]).sort_index()

        # final outlier assignment by cell-level max proba < outlier_cell_cutoff
        final_predict[final_cell_proba.values < self.outlier_cell_cutoff] = -1
        print(f'Assign {(final_predict == -1).sum()} cells as outliers '
              f'based on cell-level max prediction probability < {self.outlier_cell_cutoff:.3f}')

        self.label = final_predict.apply(lambda i: f'c{i:03d}' if i != -1 else self.outlier_label).values
        self.label_proba = final_cell_proba.values
        self.final_accuracy = final_acc
        return

    def save(self, output_path):
        joblib.dump(self, output_path)

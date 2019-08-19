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
from sklearn.feature_selection import RFECV
from sklearn.metrics import fbeta_score, make_scorer, confusion_matrix, balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold


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


class ConsensusClustering:
    def __init__(self, x):
        self.n_obs, self.n_pcs = x.shape
        self.X = x
        self.supervise_X = None

        adata = anndata.AnnData(X=None,
                                obs=pd.DataFrame([], index=[f'obs{i}' for i in range(self.n_obs)]),
                                var=pd.DataFrame([], index=[f'var{i}' for i in range(self.n_pcs)]))
        adata.obsm['X_pca'] = self.X
        # here neighbors should only use PCs
        self.neighbors = Neighbors(adata=adata)
        self.neighbors_computed = False

        # clustering
        self.leiden_result_df = None
        self.hdbscan = None
        self.consensus_clusters = None
        self.consensus_clusters_rescued = None
        self.consensus_clusters_proba = None

        # model training and outlier rescue
        self.supervise_model = None
        self.training_X = None
        self.training_label = None
        self.testing_X = None
        self.testing_label = None
        self.testing_proba = None
        self.outlier_X = None
        self.outlier_proba = None
        self.confusion_matrix = None
        self.balanced_accuracy = None
        return

    def fit_predict(self,
                    n_neighbors,
                    metric='euclidean',
                    neighbor_kwds=None,
                    leiden_repeats=200,
                    seed=1,
                    leiden_resolution=1,
                    leiden_kwds=None,
                    min_cluster_size=10,
                    min_cluster_portion=None,
                    min_samples=1,
                    epsilon='auto',
                    hdbscan_kwds=None,
                    n_jobs=1):
        # Construct KNN graph
        print('Run compute_neighbors')
        kwds = dict(n_neighbors=n_neighbors,
                    metric=metric)
        if neighbor_kwds is None:
            neighbor_kwds = kwds
        else:
            neighbor_kwds.update(kwds)
        self.compute_neighbors(**neighbor_kwds)

        # repeat Leiden clustering with different random seeds
        print('Run multi_leiden_clustering')
        kwds = dict(leiden_repeats=leiden_repeats,
                    seed=seed,
                    leiden_resolution=leiden_resolution,
                    n_jobs=n_jobs)
        if leiden_kwds is None:
            leiden_kwds = kwds
        else:
            leiden_kwds.update(kwds)
        self.multi_leiden_clustering(**leiden_kwds)

        # HDBSCAN to do consensus clustering on leiden results
        print('Run hdbscan_clustering')
        kwds = dict(min_cluster_size=min_cluster_size,
                    min_cluster_portion=min_cluster_portion,
                    min_samples=min_samples,
                    epsilon=epsilon)
        if hdbscan_kwds is None:
            hdbscan_kwds = kwds
        else:
            hdbscan_kwds.update(kwds)
        self.hdbscan_clustering(**hdbscan_kwds)

    def compute_neighbors(self,
                          n_neighbors,
                          knn=True,
                          method='umap',
                          metric='euclidean',
                          metric_kwds=None,
                          seed=1):
        metric_kwds = {} if metric_kwds is None else metric_kwds

        self.neighbors.compute_neighbors(
            n_neighbors=n_neighbors, knn=knn, n_pcs=self.n_pcs, use_rep='X_pca',
            method=method, metric=metric, metric_kwds=metric_kwds,
            random_state=seed)
        self.neighbors_computed = True
        return

    def multi_leiden_clustering(self,
                                leiden_repeats=200,
                                seed=1,
                                leiden_resolution=1,
                                partition_type=None,
                                partition_kwargs=None,
                                use_weights=True,
                                n_iterations=-1,
                                n_jobs=1):
        """Modified from scanpy"""
        if not self.neighbors_computed:
            raise ValueError('Run compute_neighbors first before multi_leiden_clustering')

        # convert neighbors to igraph
        g = self.neighbors.to_igraph()

        # generate n different seeds for each single leiden partition
        np.random.seed(seed=seed)
        random_states = np.random.choice(range(99999999), size=leiden_repeats, replace=False)
        step = max(int(leiden_repeats / n_jobs), 20)
        random_state_chunks = [random_states[i: min(i + step, leiden_repeats)]
                               for i in range(0, leiden_repeats, step)]

        results = []
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
                if leiden_resolution is not None:
                    partition_kwargs['resolution_parameter'] = leiden_resolution
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
        print(f'Repeating leiden clustering found {cluster_count.min()} - {cluster_count.max()} clusters, '
              f'mean {cluster_count.mean():.1f}, std {cluster_count.std():.2f}')
        return

    def hdbscan_clustering(self,
                           min_cluster_size=10,
                           min_cluster_portion=None,
                           min_samples=1,
                           metric='hamming',
                           cluster_selection_method='eom',
                           allow_single_cluster=True,
                           epsilon=0.2):
        if min_cluster_portion is not None:
            min_cluster_size = max(min_cluster_size, self.n_obs * min_cluster_portion)
        else:
            if min_cluster_size is None:
                raise ValueError('Either min_cluster_size or min_cluster_portion should be provided')

        runner = HDBSCAN(min_cluster_size=int(min_cluster_size),
                         min_samples=int(min_samples),
                         metric=metric,
                         cluster_selection_method=cluster_selection_method,
                         allow_single_cluster=allow_single_cluster)

        if self.leiden_result_df is None:
            raise ValueError('Run multi_leiden_clustering first before hdbscan_clustering')
        runner.fit(self.leiden_result_df)
        self.hdbscan = runner
        self.reselect_clusters(epsilon=epsilon,
                               min_cluster_size=min_cluster_size)
        return

    def reselect_clusters(self, epsilon, min_cluster_size=10):
        if epsilon == 'auto':
            # the mean cluster number of leiden clustering
            target_n_clusters = int(np.round(self.leiden_result_df.apply(lambda i: i.unique().size).min()))
            for epsilon in np.arange(0.1, 0.5, 0.02):
                self.consensus_clusters = self.hdbscan.single_linkage_tree_.get_clusters(
                    cut_distance=epsilon,
                    min_cluster_size=min_cluster_size)
                n_cluster = np.unique(self.consensus_clusters).size
                if -1 in self.consensus_clusters:
                    n_cluster -= 1
                if n_cluster <= target_n_clusters:
                    break
        else:
            self.consensus_clusters = self.hdbscan.single_linkage_tree_.get_clusters(
                cut_distance=epsilon,
                min_cluster_size=min_cluster_size)

        final_size = np.unique(self.consensus_clusters).size
        if -1 in self.consensus_clusters:
            final_size -= 1
        print(f'Final consensus clustering found {final_size} clusters')

        if self.supervise_model is not None:
            print('Consensus cluster changed, the supervised model is cleared, '
                  'redo training to fit new cluster assignment.')
        self.supervise_model = None
        self.training_X = None
        self.testing_X = None
        self.outlier_X = None
        self.confusion_matrix = None
        self.balanced_accuracy = None
        return

    def _create_model(self,
                      n_estimators=500,
                      n_splits=10,
                      fbeta=1,
                      average='weighted',
                      n_jobs=1):
        estimator = BalancedRandomForestClassifier(
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
            n_jobs=n_jobs,
            random_state=0,
            verbose=0,
            warm_start=False,
            class_weight=None)

        cv = StratifiedKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=0)

        scoring = make_scorer(
            fbeta_score,
            beta=fbeta,
            average=average)

        clf = RFECV(
            estimator,
            step=3,
            min_features_to_select=1,
            cv=cv,
            scoring=scoring,
            verbose=0,
            n_jobs=10)

        self.supervise_model = clf
        return

    def supervise_training(self,
                           x=None,
                           test_portion=0.1,
                           n_estimators=500,
                           n_splits=10,
                           fbeta=1,
                           average='weighted',
                           n_jobs=1,
                           outlier_proba_cutoff=0.5,
                           confusion_merge_cutoff=0.2,
                           seed=1):
        if self.consensus_clusters is None:
            raise ValueError('Run fit_predict first to get a clustering assignment before run supervise_training.')

        n_cluster = np.unique(self.consensus_clusters[self.consensus_clusters != -1]).size
        if n_cluster == 1:
            print('There is only one cluster except for outliers, can not train supervise model on that.')
            return
        print(f'{n_cluster} input clusters')

        if x is None:
            # use the same matrix as clustering
            self.supervise_X = self.X
        else:
            # or use the provided new matrix, could be selected raw features etc.
            raise NotImplemented('This have not been fully supported yet, '
                                 'ideally should take a full feature space and '
                                 'run feature selection only using the training sample but not testing')
            # self.supervise_X = x

        # get training, testing and outlier idx
        outlier_idx = np.where(self.consensus_clusters == -1)[0]
        valid_idx = np.where(self.consensus_clusters != -1)[0]
        n_test = max(1, int(valid_idx.size * test_portion))
        np.random.seed(seed)
        test_idx = np.random.choice(valid_idx, n_test, replace=False)
        non_training_idx = np.concatenate((test_idx, outlier_idx))
        train_idx = np.array([i for i in range(self.consensus_clusters.size)
                              if i not in non_training_idx])
        # get data
        self.training_X = self.supervise_X[train_idx, :]
        self.training_label = self.consensus_clusters[train_idx]
        self.testing_X = self.supervise_X[test_idx, :]
        self.testing_label = self.consensus_clusters[test_idx]
        self.outlier_X = self.supervise_X[outlier_idx, :]
        print(f'{len(train_idx)} training obs')
        print(f'{len(test_idx)} testing obs')
        print(f'{len(outlier_idx)} outliers')

        # create model
        print('Run RFECV')
        self._create_model(
            n_estimators=n_estimators,
            n_splits=n_splits,
            fbeta=fbeta,
            average=average,
            n_jobs=n_jobs)
        # fit model
        self.supervise_model.fit(self.training_X, self.training_label)

        print('Final testing data')
        # test model
        testing_predict_label = self.supervise_model.predict(self.testing_X)
        self.testing_proba = self.supervise_model.predict_proba(self.testing_X)
        # test performance
        self.confusion_matrix = confusion_matrix(self.testing_label, testing_predict_label)
        self.balanced_accuracy = balanced_accuracy_score(self.testing_label,
                                                         testing_predict_label)
        print(f'Overall balanced accuracy: {self.balanced_accuracy:.2f}')

        # check confusion matrix, if there are confusion pairs, merge them and redo training
        confusion_groups = {}
        normalized_confusion = self.confusion_matrix / self.confusion_matrix.sum(axis=0)
        for a, b in zip(*np.where(normalized_confusion > confusion_merge_cutoff)):
            if a == b:
                continue
            added = False
            for head, member in confusion_groups.items():
                # add pair to existing group
                if (a in member) or (b in member):
                    member.add(a)
                    member.add(b)
                    added = True
            if not added:
                confusion_groups[a] = {a, b}
        if len(confusion_groups) != 0:
            print(confusion_groups)
            print(f'Found {len(confusion_groups)} groups of clusters have > {confusion_merge_cutoff} '
                  f'confusing population in testing data. Merge each group of clusters together and redo training.')
            merge_map = {}
            for head, member in confusion_groups.items():
                for c in member:
                    merge_map[c] = head
            # merge confusing cluster groups
            self.consensus_clusters = np.array([merge_map[i] if i in merge_map else i
                                                for i in self.consensus_clusters])
            # then run supervise_training again, until no confusing group remained or only 1 cluster remained.
            self.supervise_training(
                x=x,
                test_portion=test_portion,
                n_estimators=n_estimators,
                n_splits=n_splits,
                fbeta=fbeta,
                average=average,
                n_jobs=n_jobs,
                outlier_proba_cutoff=outlier_proba_cutoff,
                # increase cutoff to fasten convergence
                confusion_merge_cutoff=confusion_merge_cutoff + 0.05)
            return

        # predict_outlier in the end
        outlier_label = self.supervise_model.predict(self.outlier_X)
        self.outlier_proba = self.supervise_model.predict_proba(self.outlier_X)
        # maximum probability of a outlier prediction
        rescued = self.outlier_proba.max(axis=1) > outlier_proba_cutoff
        _outlier_idx = outlier_idx[rescued]
        _outlier_label = outlier_label[rescued]
        _clusters = self.consensus_clusters.copy()
        _clusters[_outlier_idx] = _outlier_label
        self.consensus_clusters_rescued = _clusters
        print(f'{_outlier_idx.size} / {outlier_idx.size} outliers is rescued by prediction.')
        print(f'{outlier_idx.size - _outlier_idx.size} outliers remained.')
        return

    def save_model(self, output_path):
        joblib.dump(self.supervise_model, output_path)

import joblib
import pynndescent
import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize, OneHotEncoder
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage
from .cca import cca, lsi_cca


def top_features_idx(data, n_features):
    """
    Select top features with the highest importance in CCs

    Parameters
    ----------
    data
        data.shape = (n_cc, total_features)
    n_features
        number of features to select

    Returns
    -------
    features_idx : np.array
    """
    # data.shape = (n_cc, total_features)
    n_cc = data.shape[0]
    n_features_per_dim = n_features * 10 // n_cc
    sample_range = np.arange(n_cc)[:, None]

    # get idx of n_features_per_dim features with the highest absolute loadings
    data = np.abs(data)
    idx = np.argpartition(-data, n_features_per_dim,
                          axis=1)[:, :n_features_per_dim]
    # idx.shape = (n_cc, n_features_per_dim)

    # make sure the order of first n_features_per_dim is ordered by loadings
    idx = idx[sample_range, np.argsort(-data[sample_range, idx], axis=1)]

    for i in range(n_features // n_cc + 1, n_features_per_dim):
        features_idx = np.unique(idx[:, :i].flatten())
        if len(features_idx) > n_features:
            return features_idx
    else:
        features_idx = np.unique(idx[:, :n_features_per_dim].flatten())
        return features_idx


def find_neighbor(cc1, cc2, k, random_state=0):
    """
    find all four way of neighbors for two datasets

    Parameters
    ----------
    cc1
        cc for dataset 1
    cc2
        cc for dataset 2
    k
        number of neighbors

    Returns
    -------
    11, 12, 21, 22 neighbor matrix in shape (n_cell, k)
    """
    index = pynndescent.NNDescent(cc1,
                                  metric='euclidean',
                                  n_neighbors=k + 1,
                                  random_state=random_state)
    G11 = index.neighbor_graph[0][:, 1:k + 1]
    G21 = index.query(cc2, k=k)[0]
    index = pynndescent.NNDescent(cc2,
                                  metric='euclidean',
                                  n_neighbors=k + 1,
                                  random_state=random_state)
    G22 = index.neighbor_graph[0][:, 1:k + 1]
    G12 = index.query(cc1, k=k)[0]
    return G11, G12, G21, G22


def find_mnn(G12, G21, kanchor):
    """Calculate mutual nearest neighbor for two datasets"""
    anchor = [[i, G12[i, j]]
              for i in range(G12.shape[0])
              for j in range(kanchor)
              if (i in G21[G12[i, j], :kanchor])]
    return np.array(anchor)


def min_max(tmp, q_left=1, q_right=90):
    """normalize to q_left, q_right quantile to 0, 1, and cap extreme values"""
    tmin, tmax = np.percentile(tmp, [q_left, q_right])
    tmp = (tmp - tmin) / (tmax - tmin)
    tmp[tmp > 1] = 1
    tmp[tmp < 0] = 0
    return tmp


def filter_anchor(anchor,
                  adata_ref=None,
                  adata_qry=None,
                  high_dim_feature=None,
                  k_filter=200,
                  random_state=0):
    """
    Check if an anchor is still an anchor when only
    using the high_dim_features to construct KNN graph.
    If not, remove the anchor.
    """
    ref_data = normalize(adata_ref.X[:, high_dim_feature], axis=1)
    qry_data = normalize(adata_qry.X[:, high_dim_feature], axis=1)
    index = pynndescent.NNDescent(ref_data,
                                  metric='euclidean',
                                  n_neighbors=k_filter,
                                  random_state=random_state)
    G = index.query(qry_data, k=k_filter)[0]
    input_anchors = anchor.shape[0]
    anchor = np.array([xx for xx in anchor if (xx[0] in G[xx[1]])])
    print(f'Anchor selected with high CC feature graph: {anchor.shape[0]} / {input_anchors}')
    return anchor


def score_anchor(anchor,
                 G11,
                 G12,
                 G21,
                 G22,
                 k_score=30,
                 Gp1=None,
                 Gp2=None,
                 k_local=50):
    """
    score the anchor by the number of shared neighbors

    Parameters
    ----------
    anchor
        anchor in shape (n_anchor, 2)
    G11
    G12
    G21
    G22
    k_score
        number of neighbors to score the anchor
    Gp1
        Intra-dataset1 kNN graph
    Gp2
        Intra-dataset2 kNN graph
    k_local
        number of neighbors to calculate the local score

    Returns
    -------
    anchor with score in shape (n_anchor, 3): pd.DataFrame
    """
    tmp = [
        len(set(G11[x, :k_score]).intersection(G21[y, :k_score])) +
        len(set(G12[x, :k_score]).intersection(G22[y, :k_score]))
        for x, y in anchor
    ]
    anchor_df = pd.DataFrame(anchor, columns=['x1', 'x2'])
    anchor_df['score'] = min_max(tmp)

    if k_local:
        # if k_local is not None, then use local KNN to adjust the score
        share_nn = np.array([
            len(set(Gp1[i]).intersection(G11[i, :k_local]))
            for i in range(len(Gp1))
        ])
        tmp = [share_nn[xx] for xx in anchor_df['x1'].values]
        anchor_df['score_local1'] = min_max(tmp)

        share_nn = np.array([
            len(set(Gp2[i]).intersection(G22[i, :k_local]))
            for i in range(len(Gp2))
        ])
        tmp = [share_nn[xx] for xx in anchor_df['x2'].values]
        anchor_df['score_local2'] = min_max(tmp)

        anchor_df['score'] = anchor_df['score'] * anchor_df[
            'score_local1'] * anchor_df['score_local2']
    return anchor_df


def find_order(dist, ncell):
    """
    use dendrogram to find the order of dataset pairs
    """
    D = linkage(1 / dist, method='average')
    node_dict = {i: [i] for i in range(len(ncell))}
    alignment = []
    for xx in D[:, :2].astype(int):
        if ncell[xx[0]] < ncell[xx[1]]:
            xx = xx[::-1]
        alignment.append([node_dict[xx[0]], node_dict[xx[1]]])
        node_dict[len(ncell)] = node_dict[xx[0]] + node_dict[xx[1]]
        ncell.append(ncell[xx[0]] + ncell[xx[1]])
    return alignment


def find_nearest_anchor(data,
                        data_qry,
                        all_anchor,
                        ref,
                        qry,
                        key_correct='X_pca',
                        npc=30,
                        kweight=100,
                        sd=1,
                        random_state=0):
    print('Initialize')
    cum_ref, cum_qry = [0], [0]
    for xx in ref:
        cum_ref.append(cum_ref[-1] + data[xx].shape[0])
    for xx in qry:
        cum_qry.append(cum_qry[-1] + data[xx].shape[0])

    anchor = []
    for i, xx in enumerate(ref):
        for j, yy in enumerate(qry):
            if xx < yy:
                tmp = all_anchor[(xx, yy)].copy()
            else:
                tmp = all_anchor[(yy, xx)].copy()
                tmp[['x1', 'x2']] = tmp[['x2', 'x1']]
            tmp['x1'] += cum_ref[i]
            tmp['x2'] += cum_qry[j]
            anchor.append(tmp)
    anchor = pd.concat(anchor)
    score = anchor['score'].values
    anchor = anchor[['x1', 'x2']].values

    if key_correct == 'X':
        model = PCA(n_components=npc,
                    svd_solver='arpack',
                    random_state=random_state)
        reduce_qry = model.fit_transform(data_qry)
    else:
        reduce_qry = data_qry

    print('Find nearest anchors')
    index = pynndescent.NNDescent(reduce_qry[anchor[:, 1]],
                                  metric='euclidean',
                                  n_neighbors=kweight,
                                  random_state=random_state)
    G, D = index.query(reduce_qry, k=kweight)

    print('Normalize graph')
    cell_filter = (D[:, -1] == 0)
    D = (1 - D / D[:, -1][:, None]) * score[G]
    D[cell_filter] = score[G[cell_filter]]
    D = 1 - np.exp(-D * (sd ** 2) / 4)
    D = D / (np.sum(D, axis=1) + 1e-6)[:, None]
    return anchor, G, D, cum_qry


def transform(data,
              all_anchor,
              ref,
              qry,
              npc=30,
              k_weight=100,
              sd=1,
              chunk_size=50000,
              random_state=0,
              row_normalize=True, ):
    data_ref = np.concatenate(data[ref])
    data_qry = np.concatenate(data[qry])
    anchor, G, D, cum_qry = find_nearest_anchor(
        data=data,
        data_qry=data_qry,
        all_anchor=all_anchor,
        ref=ref,
        qry=qry,
        npc=npc,
        kweight=k_weight,
        sd=sd,
        random_state=random_state)

    print('Transform data')
    bias = data_ref[anchor[:, 0]] - data_qry[anchor[:, 1]]
    data_prj = np.zeros(data_qry.shape)

    for chunk_start in np.arange(0, data_prj.shape[0], chunk_size):
        data_prj[chunk_start:(chunk_start + chunk_size)] = \
            data_qry[chunk_start:(chunk_start + chunk_size)] + \
            (D[chunk_start:(chunk_start + chunk_size), :, None] *
             bias[G[chunk_start:(chunk_start + chunk_size)]]).sum(axis=1)
    for i, xx in enumerate(qry):
        _data = data_prj[cum_qry[i]:cum_qry[i + 1]]
        if row_normalize:
            _data = normalize(_data, axis=1)
        data[xx] = _data
    return data


class SeuratIntegration:
    def __init__(self, random_state=0):
        # intra-dataset KNN graph
        self.k_local = None
        self.key_local = None
        self.local_knn = []

        self.adata_list = []
        self.n_dataset = 0
        self.n_cells = []
        self.alignments = None
        self.all_pairs = np.array([])
        self._get_all_pairs()

        self.anchor = {}
        self.mutual_knn = {}
        self.raw_anchor = {}

        self.random_state = random_state
        pass

    def __repr__(self):
        # TODO
        pass

    def _calculate_local_knn(self):
        """If klocal is provided, we calculate the local knn graph to
        evaluate whether the anchor preserves local structure within the dataset.
        One can use a different obsm with key_local to compute knn for each dataset.
        """
        if self.k_local is not None:
            print('Find neighbors within datasets')
            for adata in self.adata_list:
                index = pynndescent.NNDescent(adata.obsm[self.key_local],
                                              metric='euclidean',
                                              n_neighbors=self.k_local + 1,
                                              random_state=self.random_state)
                self.local_knn.append(index.neighbor_graph[0][:, 1:])
        else:
            self.local_knn = [None for _ in self.adata_list]

    def _get_all_pairs(self):
        if self.alignments is not None:
            all_pairs = []
            for pair in self.alignments:
                for xx in pair[0]:
                    for yy in pair[1]:
                        if xx < yy:
                            all_pairs.append(f'{xx}-{yy}')
                        else:
                            all_pairs.append(f'{yy}-{xx}')
            self.all_pairs = np.unique(all_pairs)
        else:
            self.all_pairs = np.array([])

    def _prepare_matrix(self, i, j, key_anchor):
        adata_list = self.adata_list

        if key_anchor == 'X':
            # in case the adata var is not in the same order
            # select and order the var to make sure it is matched
            if (adata_list[i].shape[1] != adata_list[j].shape[1]) or (
                    (adata_list[i].var.index == adata_list[j].var.index).sum() <
                    adata_list[i].shape[1]
            ):
                sel_b = adata_list[i].var.index & adata_list[j].var.index
                U = adata_list[i][:, sel_b].X.copy()
                V = adata_list[j][:, sel_b].X.copy()
            else:
                U = adata_list[i].X.copy()
                V = adata_list[j].X.copy()
        else:
            U = adata_list[i].obsm[key_anchor]
            V = adata_list[j].obsm[key_anchor]

        return U, V

    def _calculate_mutual_knn_and_raw_anchors(self, i, j, U, V, k, k_anchor):
        """
        Calculate the mutual knn graph and raw anchors and
        save the results to self.mutual_knn and self.raw_anchor
        """
        G11, G12, G21, G22 = find_neighbor(U, V, k=k)
        raw_anchors = find_mnn(G12, G21, k_anchor)
        self.mutual_knn[(i, j)] = (G11, G12, G21, G22)
        self.raw_anchor[(i, j)] = raw_anchors
        return G11, G12, G21, G22, raw_anchors

    def find_anchor(self,
                    adata_list,
                    k_local=None,
                    key_local='X_pca',
                    key_anchor='X',
                    dim_red='pca',
                    scale1=False,
                    scale2=False,
                    k_filter=None,
                    n_features=200,
                    n_components=None,
                    max_cc_cells=50000,
                    k_anchor=5,
                    k_score=30,
                    alignments=None):
        self.adata_list = adata_list
        self.n_dataset = len(adata_list)
        self.n_cells = [adata.shape[0] for adata in adata_list]

        # intra-dataset KNN for scoring the anchors
        self.k_local = k_local
        self.key_local = key_local
        self._calculate_local_knn()

        # alignments and all_pairs
        self.alignments = alignments
        self._get_all_pairs()

        print('Find anchors across datasets')
        for i in range(self.n_dataset - 1):
            for j in range(i + 1, self.n_dataset):
                if (alignments is not None) and (f'{i}-{j}' not in self.all_pairs):
                    continue

                # 1. prepare input matrix for CCA
                U, V = self._prepare_matrix(i, j, key_anchor=key_anchor)

                # 2. run cca between datasets
                print('Run CCA')
                if dim_red == 'pca':
                    U, V = cca(U, V, scale1=scale1, scale2=scale2, n_components=n_components)
                elif dim_red == 'lsi':
                    U, V = lsi_cca(U, V, n_components=n_components, max_cc_cell=max_cc_cells)
                    # U.shape = (self.n_cells[i], n_components)
                    # V.shape = (self.n_cells[j], n_components)

                # 3. normalize CCV per sample/row
                U = normalize(U, axis=1)
                V = normalize(V, axis=1)

                # 4. find MNN of U and V to find anchors
                print('Find Anchors')
                _k = max([i for i in [k_anchor, k_local, k_score, 50]
                          if i is not None])
                G11, G12, G21, G22, raw_anchors = \
                    self._calculate_mutual_knn_and_raw_anchors(
                        i=i,
                        j=j,
                        U=U,
                        V=V,
                        k=_k,
                        k_anchor=k_anchor)

                # 5. filter anchors by high dimensional neighbors
                if k_filter is not None:
                    # compute ccv feature loading
                    mat = np.concatenate([U, V], axis=0).T.dot(
                        np.concatenate([adata_list[i].X, adata_list[j].X], axis=0)
                    )
                    high_dim_feature = top_features_idx(mat, n_features=n_features)

                    if self.n_cells[i] >= self.n_cells[j]:
                        raw_anchors = filter_anchor(anchor=raw_anchors,
                                                    adata_ref=adata_list[i],
                                                    adata_qry=adata_list[j],
                                                    high_dim_feature=high_dim_feature,
                                                    k_filter=k_filter,
                                                    random_state=self.random_state)
                    else:
                        raw_anchors = filter_anchor(anchor=raw_anchors[:, ::-1],
                                                    adata_ref=adata_list[j],
                                                    adata_qry=adata_list[i],
                                                    high_dim_feature=high_dim_feature,
                                                    k_filter=k_filter,
                                                    random_state=self.random_state)[:, ::-1]

                # 6. score anchors with snn and local structure preservation
                print('Score Anchors')
                anchor_df = score_anchor(anchor=raw_anchors,
                                         G11=G11,
                                         G12=G12,
                                         G21=G21,
                                         G22=G22,
                                         k_score=k_score,
                                         k_local=k_local,
                                         Gp1=self.local_knn[i],
                                         Gp2=self.local_knn[j])

                # 7. save anchors
                self.anchor[(i, j)] = anchor_df.copy()
                print(f'Identified {len(self.anchor[i, j])} anchors between datasets {i} and {j}.')
        return

    def integrate(self,
                  key_correct,
                  row_normalize=True,
                  n_components=30,
                  k_weight=100,
                  sd=1,
                  alignments=None):
        if alignments is not None:
            self.alignments = alignments

        # find order of pairwise dataset merging with hierarchical clustering
        if self.alignments is None:
            dist = []
            for i in range(self.n_dataset - 1):
                for j in range(i + 1, self.n_dataset):
                    dist.append(len(self.anchor[(i, j)]) /
                                min([self.n_cells[i], self.n_cells[j]]))
            self.alignments = find_order(np.array(dist), self.n_cells)
            print(f'Alignments: {self.alignments}')

        print('Merge datasets')
        adata_list = self.adata_list

        # initialize corrected with original data
        if key_correct == 'X':
            # correct the original feature matrix
            corrected = [adata_list[i].X.copy()
                         for i in range(self.n_dataset)]
        else:
            # correct dimensionality reduced matrix only
            corrected = [normalize(adata_list[i].obsm[key_correct], axis=1)
                         for i in range(self.n_dataset)]

        for xx in self.alignments:
            print(xx)
            corrected = transform(data=np.array(corrected),
                                  all_anchor=self.anchor,
                                  ref=xx[0],
                                  qry=xx[1],
                                  npc=n_components,
                                  k_weight=k_weight,
                                  sd=sd,
                                  random_state=self.random_state,
                                  row_normalize=row_normalize)
        return corrected

    def label_transfer(self,
                       ref,
                       qry,
                       categorical_key=None,
                       continuous_key=None,
                       key_dist='X_pca',
                       kweight=100,
                       npc=30,
                       sd=1,
                       chunk_size=50000,
                       random_state=0):
        adata_list = self.adata_list

        data_qry = np.concatenate(
            [normalize(adata_list[i].obsm[key_dist], axis=1) for i in qry])
        data_qry_index = np.concatenate([adata_list[i].obs_names for i in qry])

        anchor, G, D, cum_qry = find_nearest_anchor(data=adata_list,
                                                    all_anchor=self.anchor,
                                                    data_qry=data_qry,
                                                    ref=ref,
                                                    qry=qry,
                                                    npc=npc,
                                                    kweight=kweight,
                                                    sd=sd,
                                                    random_state=random_state)

        print('Label transfer')
        label_ref = []
        columns = []

        if categorical_key is None:
            categorical_key = []
        if continuous_key is None:
            continuous_key = []
        if len(categorical_key) == 0 and len(continuous_key) == 0:
            raise ValueError('No categorical or continuous key specified.')

        if len(categorical_key) > 0:
            tmp = pd.concat([adata_list[i].obs[categorical_key] for i in ref],
                            axis=0)
            enc = OneHotEncoder()
            label_ref.append(
                enc.fit_transform(
                    tmp[categorical_key].values.astype(np.str_)
                ).toarray()
            )
            columns += enc.categories_

        if len(continuous_key) > 0:
            tmp = pd.concat([adata_list[i].obs[continuous_key]
                             for i in ref],
                            axis=0)
            label_ref.append(tmp[continuous_key].values)
            columns += [[xx] for xx in continuous_key]

        label_ref = np.concatenate(label_ref, axis=1)
        label_qry = np.zeros((data_qry.shape[0], label_ref.shape[1]))

        bias = label_ref[anchor[:, 0]]
        for chunk_start in np.arange(0, label_qry.shape[0], chunk_size):
            label_qry[chunk_start:(chunk_start + chunk_size)] = (
                    D[chunk_start:(chunk_start + chunk_size), :, None] *
                    bias[G[chunk_start:(chunk_start + chunk_size)]]
            ).sum(axis=1)

        label_qry = pd.DataFrame(label_qry,
                                 index=data_qry_index,
                                 columns=np.concatenate(columns))
        result = {}
        for xx, yy in zip(categorical_key + continuous_key, columns):
            result[xx] = label_qry[yy]
        return result

    def save(self,
             output_path,
             save_local_knn=False,
             save_raw_anchor=False,
             save_mutual_knn=False):
        self.adata_list = []

        if not save_local_knn:
            self.local_knn = []
        if not save_raw_anchor:
            self.raw_anchor = {}
        if not save_mutual_knn:
            self.mutual_knn = {}

        joblib.dump(self, output_path)
        return

    @classmethod
    def load(cls, input_path):
        return joblib.load(input_path)

    @classmethod
    def save_transfer_results_to_adata(cls,
                                       adata,
                                       transfer_results,
                                       new_label_suffix='_transfer'):
        for key, df in transfer_results.items():
            adata.obs[key + new_label_suffix] = adata.obs[key].copy()
            adata.obs.loc[df.index, key + new_label_suffix] = df.idxmax(axis=1).values
        return

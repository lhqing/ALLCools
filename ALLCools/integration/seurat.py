import numpy as np
import pandas as pd
import pynndescent
from scipy.cluster.hierarchy import linkage
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
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


def filter_anchor(anchor,
                  adata_ref=None,
                  adata_qry=None,
                  high_dim_feature=None,
                  kfilter=200):
    """
    Check if the anchor is also a anchor in the high dimensional space of ref cells.

    Construct a new knn graph with only the high_dim_feature,
    keep anchors if query still within kfilter neighbors of reference
    """
    ref_data = normalize(adata_ref.X[:, high_dim_feature], axis=1)
    qry_data = normalize(adata_qry.X[:, high_dim_feature], axis=1)
    index = pynndescent.NNDescent(ref_data,
                                  metric='euclidean',
                                  n_neighbors=kfilter)
    G = index.query(qry_data, k=kfilter)[0]
    input_anchors = anchor.shape[0]
    anchor = np.array([xx for xx in anchor if (xx[0] in G[xx[1]])])
    print(f'Anchor selected with high CC feature graph: {anchor.shape[0]} / {input_anchors}')
    return anchor


def min_max(tmp, q_left=1, q_right=90):
    """normalize to q_left, q_right quantile to 0, 1, and cap extreme values"""
    tmin, tmax = np.percentile(tmp, [q_left, q_right])
    tmp = (tmp - tmin) / (tmax - tmin)
    tmp[tmp > 1] = 1
    tmp[tmp < 0] = 0
    return tmp


def score_anchor(anchor,
                 G11,
                 G12,
                 G21,
                 G22,
                 kscore=30,
                 Gp1=None,
                 Gp2=None,
                 klocal=50):
    """
    score the anchor by the number of shared neighbors

    Parameters
    ----------
    anchor
    G11
    G12
    G21
    G22
    kscore
    Gp1
    Gp2
    klocal

    Returns
    -------

    """
    tmp = [
        len(set(G11[x, :kscore]).intersection(G21[y, :kscore])) +
        len(set(G12[x, :kscore]).intersection(G22[y, :kscore]))
        for x, y in anchor
    ]
    # len(temp) = len(anchor)

    anchor = pd.DataFrame(anchor, columns=['x1', 'x2'])
    anchor['score'] = min_max(tmp)

    if klocal:
        share_nn = np.array([
            len(set(Gp1[i]).intersection(G11[i, :klocal]))
            for i in range(len(Gp1))
        ])
        tmp = [share_nn[xx] for xx in anchor['x1'].values]
        anchor['score_local1'] = min_max(tmp)

        share_nn = np.array([
            len(set(Gp2[i]).intersection(G22[i, :klocal]))
            for i in range(len(Gp2))
        ])
        tmp = [share_nn[xx] for xx in anchor['x2'].values]
        anchor['score_local2'] = min_max(tmp)

        anchor['score'] = anchor['score'] * anchor['score_local1'] * anchor[
            'score_local2']

    return anchor


def transform(data,
              anchor_all,
              ref,
              qry,
              key_correct,
              npc=30,
              k_weight=100,
              sd=1):
    data_ref = np.concatenate(data[ref])
    data_qry = np.concatenate(data[qry])
    cum_ref, cum_qry = [0], [0]
    for xx in ref:
        cum_ref.append(cum_ref[-1] + data[xx].shape[0])
    for xx in qry:
        cum_qry.append(cum_qry[-1] + data[xx].shape[0])
    anchor = []
    for i, xx in enumerate(ref):
        for j, yy in enumerate(qry):
            if xx < yy:
                tmp = anchor_all[(xx, yy)].copy()
            else:
                tmp = anchor_all[(yy, xx)].copy()
                tmp[['x1', 'x2']] = tmp[['x2', 'x1']]
            tmp['x1'] += cum_ref[i]
            tmp['x2'] += cum_qry[j]
            anchor.append(tmp)
    anchor = pd.concat(anchor)
    score = anchor['score'].values

    anchor = anchor[['x1', 'x2']].values
    bias = data_ref[anchor[:, 0]] - data_qry[anchor[:, 1]]
    if key_correct == 'X':
        model = PCA(n_components=npc, svd_solver='arpack')
        reduce_qry = model.fit_transform(data_qry)
    else:
        reduce_qry = data_qry
    index = pynndescent.NNDescent(reduce_qry[anchor[:, 1]],
                                  metric='euclidean',
                                  n_neighbors=50)
    G, D = index.query(reduce_qry, k=k_weight)

    cell_filter = (D[:, -1] == 0)
    D = (1 - D / D[:, -1][:, None]) * score[G]
    D[cell_filter] = score[G[cell_filter]]
    D = 1 - np.exp(-D * (sd ** 2) / 4)
    D = D / (np.sum(D, axis=1) + 1e-6)[:, None]
    data_prj = data_qry + (D[:, :, None] * bias[G]).sum(axis=1)

    for i, xx in enumerate(qry):
        data[xx] = data_prj[cum_qry[i]:cum_qry[i + 1]]
    return data


def find_neighbor(cc1, cc2, k):
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
    index = pynndescent.NNDescent(cc1, metric='euclidean', n_neighbors=k + 1)
    G11 = index.neighbor_graph[0][:, 1:k + 1]
    G21 = index.query(cc2, k=k)[0]
    index = pynndescent.NNDescent(cc2, metric='euclidean', n_neighbors=k + 1)
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


def find_order(dist, ncell):
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


def integrate(adata_list,
              key_correct='X_pca',
              key_anchor='X',
              dimred='pca',
              scale1=False,
              scale2=False,
              ncc=30,
              kanchor=5,
              kscore=30,
              kweight=100,
              klocal=None,
              key_local='X_pca',
              kfilter=None,
              n_features=200,
              npc=30,
              sd=1,
              random_state=0,
              alignments=None):
    nds = len(adata_list)
    ncell = [xx.shape[0] for xx in adata_list]

    # if score anchor by local structure preservation, compute knn for individual dataset
    # One can use a different obsm with key_local to compute knn for each dataset
    if klocal:
        print('Find neighbors within datasets')
        Gp = []
        for i in range(nds):
            index = pynndescent.NNDescent(adata_list[i].obsm[key_local],
                                          metric='euclidean',
                                          n_neighbors=klocal + 1)
            Gp.append(index.neighbor_graph[0][:, 1:])
    else:
        Gp = [None for _ in range(nds)]

    print('Find anchors across datasets')
    dist = []
    anchor = {}

    if alignments is not None:
        allpairs = []
        for pair in alignments:
            for xx in pair[0]:
                for yy in pair[1]:
                    if xx < yy:
                        allpairs.append(f'{xx}-{yy}')
                    else:
                        allpairs.append(f'{yy}-{xx}')
        allpairs = np.unique(allpairs)
    else:
        allpairs = []

    for i in range(nds - 1):
        for j in range(i + 1, nds):
            if (alignments is not None) and (f'{i}-{j}' not in allpairs):
                continue
            # run cca between datasets
            if (key_anchor == 'X') and (dimred == 'pca'):
                U, V = cca(adata_list[i].X.copy(),
                           adata_list[j].X.copy(),
                           scale1=scale1,
                           scale2=scale2,
                           n_components=ncc,
                           random_state=random_state)
            elif (key_anchor == 'X') and (dimred == 'lsi'):
                U, V = lsi_cca(adata_list[i].X.copy(),
                               adata_list[j].X.copy(),
                               n_components=ncc,
                               random_state=random_state)
            else:
                U, V = cca(adata_list[i].obsm[key_anchor].copy(),
                           adata_list[j].obsm[key_anchor].copy(),
                           scale1=scale1,
                           scale2=scale2,
                           n_components=ncc,
                           random_state=random_state)

            # compute cca feature loading
            if kfilter:
                # cc.shape = (n_cell1+n_cell2, n_cc)
                cc = np.concatenate([U, V], axis=0)
                # matrix.shape = (n_cell1+n_cell2, total_feature)
                matrix = np.concatenate([adata_list[i].X, adata_list[j].X],
                                        axis=0)
                high_dim_feature_idx = top_features_idx(cc.T.dot(matrix), n_features=n_features)
            else:
                high_dim_feature_idx = None

            # normalize ccv
            U = normalize(U, axis=1)
            V = normalize(V, axis=1)

            # fine neighbors between and within datasets
            G11, G12, G21, G22 = find_neighbor(
                U, V, k=max([kanchor, klocal, kscore, 50]))

            # find mnn as anchors
            anchor[(i, j)] = find_mnn(G12, G21, kanchor)

            # filter anchors by high dimensional neighbors
            if kfilter:
                print(f'Filter anchors')
                # use larger dataset as ref, small as query
                if ncell[i] >= ncell[j]:
                    anchor[(i, j)] = filter_anchor(anchor=anchor[(i, j)],
                                                   adata_ref=adata_list[i],
                                                   adata_qry=adata_list[j],
                                                   high_dim_feature=high_dim_feature_idx,
                                                   kfilter=kfilter)
                else:
                    # put larger dataset first, after filtering, put the order back
                    anchor[(i, j)] = filter_anchor(anchor=anchor[(i, j)][:, ::-1],
                                                   adata_ref=adata_list[j],
                                                   adata_qry=adata_list[i],
                                                   high_dim_feature=high_dim_feature_idx,
                                                   kfilter=kfilter)[:, ::-1]

            # score anchors with snn and local structure preservation
            anchor[(i, j)] = score_anchor(anchor[(i, j)],
                                          G11,
                                          G12,
                                          G21,
                                          G22,
                                          kscore=kscore,
                                          klocal=klocal,
                                          Gp1=Gp[i],
                                          Gp2=Gp[j])
            # anchor value is a dataframe, with x, y, score

            # distance between datasets
            dist.append(len(anchor[(i, j)]) / min([ncell[i], ncell[j]]))
            print(f'Identified {len(anchor[i, j])} anchors between datasets '
                  f'{i} and {j}')

    print(dist)
    print('Merge datasets')

    # find order of pairwise dataset merging with hierarchical clustering
    if alignments is None:
        alignments = find_order(np.array(dist), ncell)
    print(alignments)

    # correct batch
    if key_correct == 'X':
        corrected = [adata_list[i].X.copy() for i in range(nds)]
    else:
        corrected = [
            normalize(adata_list[i].obsm[key_correct], axis=1)
            for i in range(nds)
        ]

    for xx in alignments:
        corrected = transform(np.array(corrected),
                              anchor,
                              xx[0],
                              xx[1],
                              key_correct,
                              npc=npc,
                              k_weight=kweight,
                              sd=sd)
    return corrected

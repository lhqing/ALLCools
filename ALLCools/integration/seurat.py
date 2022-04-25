import numpy as np
import pynndescent
from sklearn.preprocessing import normalize
from sklearn.metrics import pairwise_distances
from scipy.stats import zscore
from scipy.cluster.hierarchy import fclusterdata, linkage, fcluster
from fbpca import pca
from scipy.sparse import vstack


def top_features(data, n_features):
    n_features_perdim = n_features * 10 // data.shape[0]
    sample_range = np.arange(data.shape[0])[:, None]
    data = np.abs(data)
    idx = np.argpartition(-data, n_features_perdim, axis=1)[:, :n_features_perdim]
    idx = idx[sample_range, np.argsort(-data[sample_range, idx], axis=1)]
    for i in range(n_features//data.shape[0] + 1, n_features_perdim):
        features = np.unique(idx[:, :i].flatten())
        if len(features)>n_features:
            return features
    return features
    
def filter_anchor(anchor, adata_ref=None, adata_qry=None, highdimfeature=None, kfilter=200):
    ref_data = normalize(adata_ref.X[:, highdimfeature], axis=1)
    qry_data = normalize(adata_qry.X[:, highdimfeature], axis=1)
    index = pynndescent.NNDescent(ref_data, metric='euclidean', n_neighbors=kfilter)
    G = index.query(qry_data, k=kfilter)[0]
    anchor = [xx for xx in anchor if (xx[0] in G[xx[1]])]
    return anchor


def min_max(tmp, qleft=1, qright=90):
    
    tmin, tmax = np.percentile(tmp, [qleft, qright])
    tmp = (tmp - tmin) / (tmax - tmin)
    tmp[tmp>1] = 1
    tmp[tmp<0] = 0
    return tmp

def score_anchor(anchor, G11, G12, G21, G22, kscore=30, Gp1=None, Gp2=None, klocal=50):

    tmp = [len(set(G11[x, :kscore]).intersection(G21[y, :kscore])) + 
           len(set(G12[x, :kscore]).intersection(G22[y, :kscore])) for x,y in anchor]

    anchor = pd.DataFrame(anchor, columns=['x1', 'x2'])
    anchor['score'] = min_max(tmp)

    if klocal:
        sharenn = np.array([len(set(Gp1[i]).intersection(G11[i, :klocal])) for i in range(len(Gp1))])
        tmp = [sharenn[xx] for xx in anchor['x1'].values]
        anchor['score_local1'] = min_max(tmp)

        sharenn = np.array([len(set(Gp2[i]).intersection(G22[i, :klocal])) for i in range(len(Gp2))])
        tmp = [sharenn[xx] for xx in anchor['x2'].values]
        anchor['score_local2'] = min_max(tmp)
        
        anchor['score'] = anchor['score'] * anchor['score_local1'] * anchor['score_local2']
    
    # anchor = anchor.loc[anchor['score']>0]
    return anchor

def transform(data, anchorall, ref, qry, key_correct, npc=30, kweight=100, sd=1):
    dataref = np.concatenate(data[ref])
    dataqry = np.concatenate(data[qry])
    #print(dataref.shape, dataqry.shape)
    cumref, cumqry = [0], [0]
    for xx in ref:
        cumref.append(cumref[-1] + data[xx].shape[0])
    for xx in qry:
        cumqry.append(cumqry[-1] + data[xx].shape[0])
    #cumqry = np.array(cumqry, dtype=int)
    anchor = []
    for i,xx in enumerate(ref):
        for j,yy in enumerate(qry):
            if xx<yy:
                tmp = anchorall[(xx, yy)].copy()
            else:
                tmp = anchorall[(yy, xx)].copy()
                tmp[['x1', 'x2']] = tmp[['x2', 'x1']]
            tmp['x1'] += cumref[i]
            tmp['x2'] += cumqry[j]
            anchor.append(tmp)
    anchor = pd.concat(anchor)
    score = anchor['score'].values

    anchor = anchor[['x1','x2']].values
    bias = dataref[anchor[:,0]] - dataqry[anchor[:,1]]
    if key_correct=='X':
        model = PCA(n_components=npc, svd_solver='arpack')
        reduceqry = model.fit_transform(dataqry)
    else:
        reduceqry = dataqry
    index = pynndescent.NNDescent(reduceqry[anchor[:,1]], metric='euclidean', n_neighbors=50)
    G, D = index.query(reduceqry, k=kweight)
    dataprj = np.zeros(dataqry.shape)
    cellfilter = (D[:,-1]==0)
    D = (1 - D/D[:,-1][:,None]) * score[G]
    D[cellfilter] = score[G[cellfilter]]
    D = 1 - np.exp(-D*(sd**2)/4)
    D = D / (np.sum(D, axis=1) + 1e-6)[:,None]
    dataprj = dataqry + (D[:,:,None] * bias[G]).sum(axis=1)
    for i,xx in enumerate(qry):
        data[xx] = dataprj[cumqry[i]:cumqry[i+1]]
    return data

def cca(data1, data2, scale1=True, scale2=True, n_components=50):
    if scale1:
        data1 = zscore(data1, axis=1)
    if scale2:
        data2 = zscore(data2, axis=1)
    X = data1.dot(data2.T)
    U, s, V = pca(X, n_components, True)
    return U, V.T

def lsi_cca(data1, data2, scale_factor = 100000, n_components=50):
    
    col_sum1 = data1.sum(axis=0).A1
    col_sum2 = data2.sum(axis=0).A1
    binfilter = np.logical_and(col_sum1>5, col_sum2>5)
    data1 = data1[:, binfilter]
    data2 = data2[:, binfilter]
    
    row_sum = data1.sum(axis=1).A1.astype(int)
    col_sum = col_sum1[binfilter]
    data1.data = data1.data / np.repeat(row_sum, row_sum)
    data1.data = np.log(data1.data * scale_factor + 1)
    idf1 = np.log(1 + data1.shape[0] / col_sum)
    
    row_sum = data2.sum(axis=1).A1.astype(int)
    col_sum = col_sum2[binfilter]
    data2.data = data2.data / np.repeat(row_sum, row_sum)
    data2.data = np.log(data2.data * scale_factor + 1)
    idf2 = np.log(1 + data2.shape[0] / col_sum)

    tf = data1.multiply(idf1).dot(data2.multiply(idf2).T)
    U, s, V = pca(tf, n_components, True)
    return U, V.T

def find_neighbor(cc1, cc2, k):
    index = pynndescent.NNDescent(cc1, metric='euclidean', n_neighbors=k+1)
    G11 = index.neighbor_graph[0][:, 1:k+1]
    G21 = index.query(cc2, k=k)[0]
    index = pynndescent.NNDescent(cc2, metric='euclidean', n_neighbors=k+1)
    G22 = index.neighbor_graph[0][:, 1:k+1]
    G12 = index.query(cc1, k=k)[0]
    return G11, G12, G21, G22

def find_mnn(G12, G21, kanchor):
    anchor = [[i, G12[i,j]] for i in range(G12.shape[0]) for j in range(kanchor) if (i in G21[G12[i,j], :kanchor])]
    return np.array(anchor)


def find_order(dist, ncell):
    D = linkage(1/dist, method='average')
    nodedict = {i:[i] for i in range(len(ncell))}
    alignment = []
    #print(D)
    for xx in D[:, :2].astype(int):
        if (ncell[xx[0]] < ncell[xx[1]]):
            xx = xx[::-1]
        alignment.append([nodedict[xx[0]], nodedict[xx[1]]])
        nodedict[len(ncell)] = nodedict[xx[0]] + nodedict[xx[1]]
        ncell.append(ncell[xx[0]]+ncell[xx[1]])
        #print(xx, alignment, nodedict, ncell)
    return alignment
    
def integrate(adata_list, key_correct='X_pca', key_local='X_pca', key_anchor='X', dimred='pca', 
              kanchor=5, kscore=30, kweight=100, kfilter=200, klocal=50, 
              scale1=False, scale2=False, ncc=30, n_features=200, npc=30, sd=1):
    
    nds = len(adata_list)
    ncell = [xx.shape[0] for xx in adata_list]
    
    # if need to score anchor by local structure preservation, compute knn for individual dataset
    if klocal:
        print('Find neighbors within datasets')
        Gp = []
        for i in range(nds):
            index = pynndescent.NNDescent(adata_list[i].obsm[key_local], metric='euclidean', n_neighbors=klocal+1)
            Gp.append(index.neighbor_graph[0][:, 1:])
    else:
        Gp = [None for i in range(nds)]

    print('Find anchors across datasets')
    dist = []
    anchor = {}
    for i in range(nds-1):
        for j in range(i+1, nds):
            # run cca
            if (key_anchor=='X') and (dimred=='pca'):
                U, V = cca(adata_list[i].X.copy(), adata_list[j].X.copy(), scale1=scale1, scale2=scale2, n_components=ncc)
            elif (key_anchor=='X') and (dimred=='lsi'):
                U, V = lsi_cca(adata_list[i].X.copy(), adata_list[j].X.copy(), n_components=ncc)
            else:
                U, V = cca(adata_list[i].obsm[key_anchor].copy(), adata_list[j].obsm[key_anchor].copy(), 
                           scale1=scale1, scale2=scale2, n_components=ncc)
            
            # compute ccv feature loading
            if kfilter:
                highdimfeature = top_features(np.concatenate([U, V], axis=0).T.dot(np.concatenate([adata_list[i].X, adata_list[j].X], axis=0)), 
                                              n_features=n_features)
            
            # normalize ccv
            U = normalize(U, axis=1)
            V = normalize(V, axis=1)
            
            # find mnn as anchors
            G11, G12, G21, G22 = find_neighbor(U, V, k=max([kanchor, klocal, kscore, 50]))
            anchor[(i,j)] = find_mnn(G12, G21, kanchor)
            
            # filter anchors by high dimensional neighbors
            if kfilter:
                if ncell[i]>=ncell[j]:
                    anchor[(i,j)] = filter_anchor(anchor[(i,j)], adata_list[i], adata_list[j], 
                                                  highdimfeature=highdimfeature, kfilter=kfilter)
                else:
                    anchor[(i,j)] = filter_anchor(anchor[(i,j)][:, ::-1], adata_list[j], adata_list[i], 
                                                  highdimfeature=highdimfeature, kfilter=kfilter)[:, ::-1]
            
            # score anchors with snn and local structure preservation
            anchor[(i,j)] = score_anchor(anchor[(i,j)], G11, G12, G21, G22, kscore=kscore, klocal=klocal, Gp1=Gp[i], Gp2=Gp[j])

            # distance between datasets
            dist.append(len(anchor[(i,j)]) / min([ncell[i], ncell[j]]))
            print('Identified', len(anchor[i,j]), 'anchors between datasets', i, 'and', j)
            
    print(dist)
    print('Merge datasets')
    
    # find order of pairwise dataset merging with hierarchical clustering
    alignments = find_order(np.array(dist), ncell)
    print(alignments)
    
    # correct batch
    corrected = adata_list.copy()
    if key_correct=='X':
        corrected = [adata_list[i].X.copy() for i in range(nds)]
    else:
        corrected = [normalize(adata_list[i].obsm[key_correct], axis=1) for i in range(nds)]
    for xx in alignments:
        corrected = transform(np.array(corrected), anchor, xx[0], xx[1], key_correct, npc=npc, kweight=kweight, sd=sd)
    return corrected





 

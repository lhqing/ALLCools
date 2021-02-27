import logging
from typing import Tuple, Any

import numpy as np
from numba import jit
from pynndescent import NNDescent
from scipy import sparse
from sklearn.neighbors import NearestNeighbors as _NearestNeighbors


@jit(["float32(float64[:], float64[:])", "float32(float32[:], float32[:])"], nopython=True, cache=True)
def jensen_shannon_divergence(pk: np.ndarray, qk: np.ndarray) -> float:
    N = pk.shape[0]
    # pk = pk / np.sum(pk)
    # qk = qk / np.sum(qk)
    m = (pk + qk) / 2

    vec = np.zeros(N)
    for i in range(N):
        if pk[i] > 0 and m[i] > 0:
            vec[i] = pk[i] * np.log(pk[i] / m[i])
        elif pk[i] == 0 and m[i] >= 0:
            vec[i] = 0
        else:
            vec[i] = np.inf
    Dpm = np.sum(vec) / np.log(2)

    vec = np.zeros(N)
    for i in range(N):
        if qk[i] > 0 and m[i] > 0:
            vec[i] = qk[i] * np.log(qk[i] / m[i])
        elif qk[i] == 0 and m[i] >= 0:
            vec[i] = 0
        else:
            vec[i] = np.inf
    Dqm = np.sum(vec) / np.log(2)

    return (Dpm + Dqm) / 2


@jit(["float32(float64[:], float64[:])", "float32(float32[:], float32[:])"], nopython=True, cache=True)
def jensen_shannon_distance(pk: np.ndarray, qk: np.ndarray) -> float:
    """
    Remarks:
        pk and qk must already be normalized so that np.sum(pk) == np.sum(qk) == 1
    """
    return np.sqrt(jensen_shannon_divergence(pk, qk))


@jit(nopython=True)
def balance_knn_loop(dsi: np.ndarray,
                     dist: np.ndarray,
                     lsi: np.ndarray,
                     maxl: int,
                     k: int) -> Tuple:
    """Fast and greedy algorythm to balance a K-NN graph so that no node is the NN to more than maxl other nodes
        Arguments
        ---------
        dsi : np.ndarray  (samples, K)
            distance sorted indexes (as returned by sklearn NN)
        dist : np.ndarray  (samples, K)
            the actual distance corresponding to the sorted indexes
        lsi : np.ndarray (samples,)
            degree of connectivity (l) sorted indexes
        maxl : int
            max degree of connectivity (from others to self) allowed
        k : int
            number of neighbours in the final graph
        return_distance : bool
            wether to return distance
        Returns
        -------
        dsi_new : np.ndarray (samples, k+1)
            indexes of the NN, first column is the sample itself
        dist_new : np.ndarray (samples, k+1)
            distances to the NN
        l: np.ndarray (samples)
            l[i] is the number of connections from other samples to the sample i
    """
    assert dsi.shape[1] >= k, "sight needs to be bigger than k"
    # numba signature "Tuple((int64[:,:], float32[:, :], int64[:]))(int64[:,:], int64[:], int64, int64, bool)"
    dsi_new = -1 * np.ones((dsi.shape[0], k + 1), np.int64)  # maybe d.shape[0]
    l = np.zeros(dsi.shape[0], np.int64)
    dist_new = np.zeros(dsi_new.shape, np.float64)
    for i in range(dsi.shape[0]):  # For every node
        el = lsi[i]  # start from high degree of connectivity
        p = 0
        j = 0
        for j in range(dsi.shape[1]):  # For every other node it is connected (sight)
            if p >= k:
                break
            m = dsi[el, j]
            if el == m:
                dsi_new[el, 0] = el
                continue
            if l[m] >= maxl:
                # skip the node if it is already connected with maxl nodes
                continue
            dsi_new[el, p + 1] = m
            l[m] = l[m] + 1
            dist_new[el, p + 1] = dist[el, j]
            p += 1
        if (j == dsi.shape[1] - 1) and (p < k):
            while p < k:
                dsi_new[el, p + 1] = el
                dist_new[el, p + 1] = dist[el, 0]
                p += 1
    return dist_new, dsi_new, l


def knn_balance(dsi: np.ndarray,
                dist: np.ndarray = None,
                maxl: int = 200,
                k: int = 60) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Balance a K-NN graph so that no node is the NN to more than maxl other nodes
        Arguments
        ---------
        dsi : np.ndarray  (samples, K)
            distance sorted indexes (as returned by sklearn NN)
        dist : np.ndarray  (samples, K)
            the actual distance corresponding to the sorted indexes
        maxl : int
            max degree of connectivity allowed
        k : int
            number of neighbours in the final graph
        Returns
        -------
        dist_new : np.ndarray (samples, k+1)
            distances to the NN
        dsi_new : np.ndarray (samples, k+1)
            indexes of the NN, first column is the sample itself
        l: np.ndarray (samples)
            l[i] is the number of connections from other samples to the sample i
    """
    l = np.bincount(dsi.flat[:], minlength=dsi.shape[0])  # degree of connectivity
    lsi = np.argsort(l, kind="mergesort")[::-1]  # element index sorted by descending degree of connectivity
    return balance_knn_loop(dsi, dist, lsi, maxl, k)


class NearestNeighbors:
    """Greedy algorithm to balance a K-nearest neighbour graph
    It has an API similar to scikit-learn
    Parameters
    ----------
    k : int  (default=50)
        the number of neighbours in the final graph
    sight_k : int  (default=100)
        the number of neighbours in the initialization graph
        It correspondent to the farthest neighbour that a sample is allowed to connect to
        when no closest neighbours are allowed. If sight_k is reached then the matrix is filled
        with the sample itself
    maxl : int  (default=200)
        max degree of connectivity allowed. Avoids the presence of hubs in the graph, it is the
        maximum number of neighbours that are allowed to contact a node before the node is blocked
    mode : str (default="connectivity")
        decide wich kind of utput "distance" or "connectivity"
    n_jobs : int  (default=4)
        parallelization of the standard KNN search preformed at initialization
    """

    def __init__(self,
                 k: int = 50,
                 sight_k: int = 100,
                 maxl: int = 200,
                 mode: str = "distance",
                 metric: str = "euclidean",
                 minkowski_p: int = 20,
                 n_jobs: int = -1) -> None:
        # input parameters
        self.k = k
        self.sight_k = sight_k
        self.maxl = maxl
        self.mode = mode
        self.metric = metric
        self.minkowski_p = minkowski_p
        self.n_jobs = n_jobs

        # NN graphs
        self.data = None
        self._nn = None  # raw KNN
        self.bknn = None  # balanced KNN
        self.dist = None  # balanced KNN distances
        self.dsi = None  # balanced KNN neighbor index
        self.l = None  # balanced KNN degree of connectivity
        self.mknn = None  # mutual KNN based on bknn
        self.rnn = None  # radius NN based on mknn

    @property
    def n_samples(self) -> int:
        return self.data.shape[0]

    def fit(self, data: np.ndarray, sight_k: int = None) -> Any:
        """Fits the model
        data: np.ndarray (samples, features)
            np
        sight_k: int
            the farthest point that a node is allowed to connect to when its closest neighbours are not allowed
        """
        self.data = data
        if sight_k is not None:
            self.sight_k = sight_k
        logging.debug(f"First search the {self.sight_k} nearest neighbours for {self.n_samples}")
        np.random.seed(13)
        if self.metric == "correlation":
            self._nn = _NearestNeighbors(n_neighbors=self.sight_k + 1, metric=self.metric, p=self.minkowski_p,
                                         n_jobs=self.n_jobs, algorithm="brute")
            self._nn.fit(self.data)
        elif self.metric == "js":
            self._nn = NNDescent(data=self.data, metric=jensen_shannon_distance)
        else:
            self._nn = _NearestNeighbors(n_neighbors=self.sight_k + 1, metric=self.metric, p=self.minkowski_p,
                                         n_jobs=self.n_jobs, leaf_size=30)
            self._nn.fit(self.data)

        # call this to calculate bknn
        self.kneighbors_graph(mode='distance')
        return self

    def kneighbors(self, X: np.ndarray = None, maxl: int = None, mode: str = "distance") -> Tuple[
        np.ndarray, np.ndarray, np.ndarray]:
        if self._nn is None:
            raise ValueError('must fit() before generating kneighbors graphs')
        """Finds the K-neighbors of a point.
            Returns indices of and distances to the neighbors of each point.
            Parameters
            ----------
            X : array-like, shape (n_query, n_features),
                The query point or points.
                If not provided, neighbors of each indexed point are returned.
                In this case, the query point is not considered its own neighbor.
            maxl: int
                max degree of connectivity allowed
            mode : "distance" or "connectivity"
                Decides the kind of output
            Returns
            -------
            dist_new : np.ndarray (samples, k+1)
                distances to the NN
            dsi_new : np.ndarray (samples, k+1)
                indexes of the NN, first column is the sample itself
            l: np.ndarray (samples)
                l[i] is the number of connections from other samples to the sample i
            NOTE:
            First column (0) correspond to the sample itself, the nearest neighbour is at the second column (1)
        """
        if X is not None:
            self.data = X
        if maxl is not None:
            self.maxl = maxl
        if mode == "distance":
            if self.metric == "js":
                self.dsi, self.dist = self._nn.query(self.data, k=self.sight_k + 1)
            else:
                self.dist, self.dsi = self._nn.kneighbors(self.data, return_distance=True)
        else:
            if self.metric == "js":
                self.dsi, _ = self._nn.query(self.data, k=self.sight_k + 1)
            else:
                self.dsi = self._nn.kneighbors(self.data, return_distance=False)
            self.dist = np.ones_like(self.dsi, dtype='float64')
            self.dist[:, 0] = 0
        logging.debug(
            f"Using the initialization network to find a {self.k}-NN "
            f"graph with maximum connectivity of {self.maxl}")
        self.dist, self.dsi, self.l = knn_balance(self.dsi, self.dist, maxl=self.maxl, k=self.k)
        return self.dist, self.dsi, self.l

    def kneighbors_graph(self, X: np.ndarray = None, maxl: int = None, mode: str = "distance") -> sparse.csr_matrix:
        """Retrun the K-neighbors graph as a sparse csr matrix
            Parameters
            ----------
            X : array-like, shape (n_query, n_features),
                The query point or points.
                If not provided, neighbors of each indexed point are returned.
                In this case, the query point is not considered its own neighbor.
            maxl: int
                max degree of connectivity allowed
            mode : "distance" or "connectivity"
                Decides the kind of output
            Returns
            -------
            neighbor_graph : scipy.sparse.csr_matrix
                The values are either distances or connectivity dependig of the mode parameter
            NOTE: The diagonal will be zero even though the value 0 is actually stored
        """
        dist_new, dsi_new, _ = self.kneighbors(X=X, maxl=maxl, mode=mode)
        logging.debug("Returning sparse matrix")
        self.bknn = sparse.csr_matrix((np.ravel(dist_new), np.ravel(dsi_new),
                                       np.arange(0, dist_new.shape[0] * dist_new.shape[1] + 1, dist_new.shape[1])),
                                      (self.n_samples, self.n_samples))
        self.bknn.eliminate_zeros()
        return self.bknn

    def mnn_graph(self):
        """get mutual nearest neighbor graph from bknn"""
        if self.mknn is None:
            if self.bknn is None:
                raise ValueError('must fit() before generating kneighbors graphs')
            # element-wise minimum between bknn and bknn.T, so non-mutual value will be 0
            self.mknn = self.bknn.minimum(self.bknn.transpose())
        return self.mknn

    def rnn_graph(self):
        """get rnn from mknn, return a sparse binary matrix"""
        # Convert distances to similarities
        if self.mknn is None:
            self.mnn_graph()
        mknn_sim = self.mknn.copy()
        bknn_sim = self.bknn.copy()
        max_d = self.bknn.data.max()
        bknn_sim.data = (max_d - bknn_sim.data) / max_d
        mknn_sim.data = (max_d - mknn_sim.data) / max_d
        mknn_sim = mknn_sim.tocoo()
        mknn_sim.setdiag(0)

        # Compute the effective resolution
        d = 1 - bknn_sim.data
        radius = np.percentile(d, 90)
        logging.info(f"  90th percentile radius: {radius:.02}")
        inside = mknn_sim.data > 1 - radius
        self.rnn = sparse.coo_matrix((mknn_sim.data[inside], (mknn_sim.row[inside], mknn_sim.col[inside])),
                                     shape=mknn_sim.shape)
        return self.rnn

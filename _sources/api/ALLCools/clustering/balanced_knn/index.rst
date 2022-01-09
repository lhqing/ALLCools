:py:mod:`ALLCools.clustering.balanced_knn`
==========================================

.. py:module:: ALLCools.clustering.balanced_knn


Module Contents
---------------

.. py:function:: jensen_shannon_divergence(pk: numpy.ndarray, qk: numpy.ndarray) -> float


.. py:function:: jensen_shannon_distance(pk: numpy.ndarray, qk: numpy.ndarray) -> float

   Remarks:
       pk and qk must already be normalized so that np.sum(pk) == np.sum(qk) == 1


.. py:function:: balance_knn_loop(dsi: numpy.ndarray, dist: numpy.ndarray, lsi: numpy.ndarray, maxl: int, k: int) -> Tuple

   Fast and greedy algorythm to balance a K-NN graph so that no node is the NN to more than maxl other nodes
   :param dsi: distance sorted indexes (as returned by sklearn NN)
   :type dsi: np.ndarray  (samples, K)
   :param dist: the actual distance corresponding to the sorted indexes
   :type dist: np.ndarray  (samples, K)
   :param lsi: degree of connectivity (l) sorted indexes
   :type lsi: np.ndarray (samples,)
   :param maxl: max degree of connectivity (from others to self) allowed
   :type maxl: int
   :param k: number of neighbours in the final graph
   :type k: int
   :param return_distance: wether to return distance
   :type return_distance: bool

   :returns: * **dsi_new** (*np.ndarray (samples, k+1)*) -- indexes of the NN, first column is the sample itself
             * **dist_new** (*np.ndarray (samples, k+1)*) -- distances to the NN
             * **l** (*np.ndarray (samples)*) -- l[i] is the number of connections from other samples to the sample i


.. py:function:: knn_balance(dsi: numpy.ndarray, dist: numpy.ndarray = None, maxl: int = 200, k: int = 60) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]

   Balance a K-NN graph so that no node is the NN to more than maxl other nodes
   :param dsi: distance sorted indexes (as returned by sklearn NN)
   :type dsi: np.ndarray  (samples, K)
   :param dist: the actual distance corresponding to the sorted indexes
   :type dist: np.ndarray  (samples, K)
   :param maxl: max degree of connectivity allowed
   :type maxl: int
   :param k: number of neighbours in the final graph
   :type k: int

   :returns: * **dist_new** (*np.ndarray (samples, k+1)*) -- distances to the NN
             * **dsi_new** (*np.ndarray (samples, k+1)*) -- indexes of the NN, first column is the sample itself
             * **l** (*np.ndarray (samples)*) -- l[i] is the number of connections from other samples to the sample i


.. py:class:: NearestNeighbors(k: int = 50, sight_k: int = 100, maxl: int = 200, mode: str = 'distance', metric: str = 'euclidean', minkowski_p: int = 20, n_jobs: int = -1)

   Greedy algorithm to balance a K-nearest neighbour graph
   It has an API similar to scikit-learn
   :param k: the number of neighbours in the final graph
   :type k: int  (default=50)
   :param sight_k: the number of neighbours in the initialization graph
                   It correspondent to the farthest neighbour that a sample is allowed to connect to
                   when no closest neighbours are allowed. If sight_k is reached then the matrix is filled
                   with the sample itself
   :type sight_k: int  (default=100)
   :param maxl: max degree of connectivity allowed. Avoids the presence of hubs in the graph, it is the
                maximum number of neighbours that are allowed to contact a node before the node is blocked
   :type maxl: int  (default=200)
   :param mode: decide wich kind of utput "distance" or "connectivity"
   :type mode: str (default="connectivity")
   :param n_jobs: parallelization of the standard KNN search preformed at initialization
   :type n_jobs: int  (default=4)

   .. py:method:: n_samples(self) -> int
      :property:


   .. py:method:: fit(self, data: numpy.ndarray, sight_k: int = None) -> Any

      Fits the model
      data: np.ndarray (samples, features)
          np
      sight_k: int
          the farthest point that a node is allowed to connect to when its closest neighbours are not allowed


   .. py:method:: kneighbors(self, X: numpy.ndarray = None, maxl: int = None, mode: str = 'distance') -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]


   .. py:method:: kneighbors_graph(self, X: numpy.ndarray = None, maxl: int = None, mode: str = 'distance') -> scipy.sparse.csr_matrix

      Retrun the K-neighbors graph as a sparse csr matrix
      :param X: The query point or points.
                If not provided, neighbors of each indexed point are returned.
                In this case, the query point is not considered its own neighbor.
      :type X: array-like, shape (n_query, n_features),
      :param maxl: max degree of connectivity allowed
      :type maxl: int
      :param mode: Decides the kind of output
      :type mode: "distance" or "connectivity"

      :returns: * **neighbor_graph** (*scipy.sparse.csr_matrix*) -- The values are either distances or connectivity dependig of the mode parameter
                * **NOTE** (*The diagonal will be zero even though the value 0 is actually stored*)


   .. py:method:: mnn_graph(self)

      get mutual nearest neighbor graph from bknn


   .. py:method:: rnn_graph(self)

      get rnn from mknn, return a sparse binary matrix




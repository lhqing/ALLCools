:py:mod:`ALLCools.pseudo_cell.pseudo_cell_knn`
==============================================

.. py:module:: ALLCools.pseudo_cell.pseudo_cell_knn


Module Contents
---------------

.. py:class:: SamplingStrategyEnum

   Bases: :py:obj:`enum.Enum`

   Generic enumeration.

   Derive from this class to define new enumerations.

   .. py:attribute:: SPARSE
      :annotation: :str = sparse

      

   .. py:attribute:: DENSE
      :annotation: :str = dense

      


.. py:class:: ExamplerAndNeighborhoodSampler(data, n_components=30, normalize=False)

   .. py:method:: _sample_pulp_dist(n_kernels, pulp_size)


   .. py:method:: _select_dist_thresh(pulp_size, n_tests=100, pulp_thicken_ratio=1.2, robust_quantile=0.9)

      Select a distance threshold to be used in the sampling procedure.

      :param pulp_size: Size of the neighborhood
      :type pulp_size: int
      :param n_tests: Sampling number
      :type n_tests: int
      :param pulp_thicken_ratio: to get a more robust distance estimate, the neighborhood size
                                 will be relaxed by this ratio during sampling
      :type pulp_thicken_ratio: float
      :param robust_quantile: to get a more robust distance estimate, this quantile of the all
                              sampled distance will be used as the final distance
      :type robust_quantile: float

      :returns: **dist_thresh** -- a general distance threshold for neighborhood of give size by random sampling.
      :rtype: str


   .. py:method:: _sample_fruit(pulp_size: int, n_kernels: int = None, max_attempts: int = 100, dist_thresh: float = None, ovlp_tol: float = 0.2, min_pulp_size: int = None, k=1000, strategy: SamplingStrategyEnum = 'dense')

      Return a list of examplers and their neighborhoods.

      Iteratively generating new examplers and corresponding neighborhoods.
      In each iteration, at most one new exampler and the neighborhood is added to the final list.

              Parameters:
                      k (int):
                      strategy (['sparse','dense']):
              Returns:
                      kernels (list): indexes of selected examplers
                      pulps (list of list): list of indexes of the corresponding neighborhoods

      :param pulp_size: desired size of the neighborhood
      :type pulp_size: int
      :param n_kernels: desired number of examplers; as many as posible if None
      :type n_kernels: int
      :param max_attempts: max number of iterations of find a new exampler
      :type max_attempts: int
      :param dist_thresh: min distance allowed between examplers
      :type dist_thresh: float
      :param ovlp_tol: max overlapping ratio between examplers
      :type ovlp_tol: float
      :param min_pulp_size: min size of the neighborhood.
                            For some examplers, the neighborhood sizes may not reach the desired size.
                            This is the mininum acceptable neighobrhood size.
      :type min_pulp_size: int
      :param k: Numbers of new exampler candidates to consider in each iteration.
                No need to change this argument in general.
      :type k: int
      :param strategy: Strategy of adding a new exampler. Two strategies are available:
                       1. 'dense' will choose the candidate closest to the examplers already-selected.
                       Recomand to combine 'dense' strategy with 'n_kernels = None'
                       2. 'sparse' will choose the candidate fartherest to the examplers already-selected.
                       Recomand to use 'sparse' strategy when a specific 'n_kernels' is desired.
      :type strategy: str


   .. py:method:: sample_examplers_and_neighborhoods(n_examplers: int, n_neighbors: int, min_n_neighbors: int = None, ovlp_tol: float = 0, dist_thresh: float = None, strategy: SamplingStrategyEnum = 'dense', max_attempts: int = 100)

      Sample examplers and their neighborhoods.



.. py:function:: sample_pseudo_cells(cell_meta, cluster_col, coords, target_pseudo_size, min_pseudo_size=None, ignore_small_cluster=False, n_components=30, pseudo_ovlp=0, n_pseudos=None, strategy: SamplingStrategyEnum = 'dense')

   Sample pseudo cells.


.. py:function:: generate_pseudo_cells(adata, cluster_col='leiden', obsm='X_pca', target_pseudo_size=100, min_pseudo_size=None, ignore_small_cluster=False, n_components=None, aggregate_func='downsample', pseudo_ovlp=0)



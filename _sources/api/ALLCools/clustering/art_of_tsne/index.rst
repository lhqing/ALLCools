:py:mod:`ALLCools.clustering.art_of_tsne`
=========================================

.. py:module:: ALLCools.clustering.art_of_tsne


Module Contents
---------------

.. py:function:: art_of_tsne(X: numpy.ndarray, metric: Union[str, Callable] = 'euclidean', exaggeration: float = -1, perplexity: int = 30, n_jobs: int = -1) -> openTSNE.TSNEEmbedding

   Implementation of Dmitry Kobak and Philipp Berens
   "The art of using t-SNE for single-cell transcriptomics" based on openTSNE.
   See https://doi.org/10.1038/s41467-019-13056-x | www.nature.com/naturecommunications
   :param X                               The data matrix of shape:
   :type X                               The data matrix of shape: n_cells, n_genes) i.e. (n_samples, n_features
   :param metric                  Any metric allowed by PyNNDescent (default: 'euclidean')
   :param exaggeration    The exaggeration to use for the embedding:
   :param perplexity              The perplexity to use for the embedding:

   :returns: The embedding as an opentsne.TSNEEmbedding object (which can be cast to an np.ndarray)


.. py:function:: tsne(adata, obsm='X_pca', metric: Union[str, Callable] = 'euclidean', exaggeration: float = -1, perplexity: int = 30, n_jobs: int = -1)



:py:mod:`ALLCools.clustering.art_of_tsne`
=========================================

.. py:module:: ALLCools.clustering.art_of_tsne

.. autoapi-nested-parse::

   The function art_of_tsne is from cytograph2 package.

   https://github.com/linnarsson-lab/cytograph2/blob/master/cytograph/embedding/art_of_tsne.py

   The idea behind that is based on :cite:p:`Kobak2019` with T-SNE algorithm implemented in
   [openTSNE](https://opentsne.readthedocs.io/en/latest/) :cite:p:`Policar2019`.



Module Contents
---------------

.. py:function:: art_of_tsne(X: numpy.ndarray, metric: Union[str, Callable] = 'euclidean', exaggeration: float = -1, perplexity: int = 30, n_jobs: int = -1) -> openTSNE.TSNEEmbedding

   Calculate T-SNE embedding with the openTSNE package.

   Implementation of Dmitry Kobak and Philipp Berens
   "The art of using t-SNE for single-cell transcriptomics" based on openTSNE.
   See https://doi.org/10.1038/s41467-019-13056-x | www.nature.com/naturecommunications

   :param X: The data matrix of shape (n_cells, n_genes) i.e. (n_samples, n_features)
   :param metric: Any metric allowed by PyNNDescent (default: 'euclidean')
   :param exaggeration: The exaggeration to use for the embedding
   :param perplexity: The perplexity to use for the embedding
   :param n_jobs: Number of CPUs to use

   :rtype: The embedding as an opentsne.TSNEEmbedding object (which can be cast to an np.ndarray)


.. py:function:: tsne(adata, obsm='X_pca', metric: Union[str, Callable] = 'euclidean', exaggeration: float = -1, perplexity: int = 30, n_jobs: int = -1)

   Calculate T-SNE embedding with the openTSNE package.

   Use the openTSNE package :cite:p:`Policar2019`
   and parameter optimization strategy described in :cite:p:`Kobak2019`.

   :param adata: adata object with principle components or equivalent matrix stored in .obsm
   :param obsm: name of the matrix in .obsm that can be used as T-SNE input
   :param metric: Any metric allowed by PyNNDescent (default: 'euclidean')
   :param exaggeration: The exaggeration to use for the embedding
   :param perplexity: The perplexity to use for the embedding
   :param n_jobs: Number of CPUs to use

   :rtype: T-SNE embedding will be stored at adata.obsm["X_tsne"]



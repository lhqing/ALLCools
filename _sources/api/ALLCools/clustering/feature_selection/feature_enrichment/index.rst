:py:mod:`ALLCools.clustering.feature_selection.feature_enrichment`
==================================================================

.. py:module:: ALLCools.clustering.feature_selection.feature_enrichment


Module Contents
---------------

.. py:function:: _calculate_enrichment_score(raw_adata, labels)

   Enrichment score modified from :cite:p:`Zeisel2018` for normalized methylation fractions


.. py:function:: _calculate_enrichment_score_cytograph(adata, labels)

   The original CEF algorithm from :cite:p:`Zeisel2018` for count based data (RNA, ATAC)


.. py:function:: _plot_enrichment_result(qvals, enrichment, null_enrichment, alpha)

   Make some plots for the p-value, q-value and # of CEFs distribution


.. py:function:: _aggregate_enrichment(adata, enrichment, top_n, alpha, qvals, cluster_col)

   Aggregate enrichment results, calculate q values


.. py:function:: cluster_enriched_features(adata: anndata.AnnData, cluster_col: str, top_n=200, alpha=0.05, stat_plot=True, method='mc')

   Calculate top Cluster Enriched Features (CEF) from per-cell normalized dataset.
   An post-clustering feature selection step adapted from :cite:p:`Zeisel2018,La_Manno2021`
   and their great [cytograph2](https://github.com/linnarsson-lab/cytograph2) package.
   For details about CEF calculation, read the methods of :cite:p:`Zeisel2018`. Note that
   in original paper, they look for cluster-specific highly expressed genes as CEFs;
   for methylation, we are looking for hypo-methylation as CEFs, so the score and test is reversed.

   :param adata: adata containing per-cell normalized values.
                 For methylation fraction, the value need to be 1-centered
                 (1 means cell's average methylation), like those produced by
                 :func:`ALLCools.mcds.mcds.MCDS.add_mc_frac` with `normalize_per_cell=True`.
                 For RNA and ATAC, you can use per cell normalized counts.
                 Do not log transform the data before running this function
   :param cluster_col: The name of categorical variable in adata.obs
   :param top_n: Select top N CEFs for each cluster
   :param alpha: FDR corrected q-value cutoff
   :param stat_plot: Whether making some summary plots for the CEF calculation
   :param method: "mc" for methylation CEF (look for hypo-methylation),
                  "rna" and "atac" for the RNA and ATAC or any count based data
                  (use the original cytograph algorithm, look for higher value)

   :returns: * *Modify adata inplace, adding a dictionary in adata.uns called f"{cluster_col}_feature_enrichment"*
             * *The dictionary contains "qvals" (np.ndarray cluster-by-feature enrichment score q-value) and*
             * *"cluster_order" (cluster order of the "qvals")*



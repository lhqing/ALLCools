:py:mod:`ALLCools.clustering.feature_selection.feature_enrichment`
==================================================================

.. py:module:: ALLCools.clustering.feature_selection.feature_enrichment


Module Contents
---------------

.. py:function:: calculate_enrichment_score(raw_adata, labels)

   Enrichment score modified from Zeisel et al. 2018 Cell for normalized methylation fractions


.. py:function:: calculate_enrichment_score_cytograph(adata, labels)

   Enrichment score algorithm from Zeisel et al. 2018 Cell


.. py:function:: _plot_enrichment_result(qvals, enrichment, null_enrichment, alpha)


.. py:function:: _aggregate_enrichment(adata, enrichment, top_n, alpha, qvals, cluster_col)


.. py:function:: cluster_enriched_features(adata, cluster_col, top_n=200, alpha=0.05, stat_plot=True, method='mc')

   :param adata:
   :param cluster_col:
   :param top_n:
   :param alpha:
   :param stat_plot:
   :param method: "mc" for methylation fraction,
                  "rna" and "atac" for the original algorithm for the RNA



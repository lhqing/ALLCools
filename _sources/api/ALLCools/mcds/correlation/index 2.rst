:py:mod:`ALLCools.mcds.correlation`
===================================

.. py:module:: ALLCools.mcds.correlation

.. autoapi-nested-parse::

   See here https://stackoverflow.com/questions/52371329/fast-spearman-correlation-between-two-pandas-dataframes

   Calculate correlation between two matrix, row by row



Module Contents
---------------

.. py:function:: _mean(a)


.. py:function:: _std(a)


.. py:function:: _corr(a, b, row, col)

   Correlation between rows in a and b, no nan value


.. py:function:: _corr_preprocess(da, sample_mch, sample_mcg, cpu=1)


.. py:function:: corr(adata_a, adata_b, max_dist, cpu=1, method='pearson', shuffle_sample=None, calculate_n=None)


.. py:function:: _verify_correlation_da(data: xarray.DataArray, region_dim, pos)


.. py:function:: region_correlation(data_a, data_b, sample_mch, sample_mcg, region_dim_a=None, region_dim_b=None, pos_a=None, pos_b=None, method='pearson', max_dist=1000000, cpu=1, null='sample', null_n=100000, chroms=None)


.. py:function:: _select_corr_cutoff(true: numpy.ndarray, null: numpy.ndarray, alpha=0.05, direction='+')


.. py:function:: get_corr_table(total_results, null_results, region_dim_a, region_dim_b, direction='-', alpha=0.05)



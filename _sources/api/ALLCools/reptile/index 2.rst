:py:mod:`ALLCools.reptile`
==========================

.. py:module:: ALLCools.reptile


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   reptile/index.rst


Package Contents
----------------

.. py:class:: REPTILE(output_path, train_regions, dmr_regions, train_region_labels, train_sample, bigwig_table, chrom_size_path, window_size=2000, step_size=200, dmr_slop=150, fillna_by_zero=None)

   .. py:method:: generate_region_ds(self)


   .. py:method:: train_region_ds(self)
      :property:


   .. py:method:: train_dmr_ds(self)
      :property:


   .. py:method:: query_region_ds(self)
      :property:


   .. py:method:: query_dmr_ds(self)
      :property:


   .. py:method:: region_model(self)
      :property:


   .. py:method:: dmr_model(self)
      :property:


   .. py:method:: _validate_region_name(self, name)


   .. py:method:: annotate_by_bigwigs(self, name, slop, cpu, redo=False)


   .. py:method:: _filter_na_train(self, name, sample, max_na_rate=0.5)


   .. py:method:: prepare_training_input(self, name)


   .. py:method:: auto_ml(data, label, output_path, train_size=0.75, random_state=42, cpu=1, tpot_generations=5, tpot_max_time_mins=60, **tpot_kwargs)
      :staticmethod:


   .. py:method:: _train(self, region_dim, slop, cpu, **kwargs)


   .. py:method:: train_region_model(self, slop=None, cpu=1, **kwargs)


   .. py:method:: train_dmr_model(self, slop=None, cpu=1, **kwargs)


   .. py:method:: fit(self, cpu=10, **kwargs)

      Convenient function to train everything by default parameters


   .. py:method:: _predict(self, region_dim, cpu, mask_cutoff)


   .. py:method:: predict(self, cpu, mask_cutoff=0.3, bw_bin_size=10, enhancer_cutoff=0.7)


   .. py:method:: _dump_sample(self, sample, mask_cutoff, bw_bin_size)


   .. py:method:: dump_bigwigs(self, cpu, mask_cutoff, bw_bin_size)


   .. py:method:: call_enhancers(self, threshold=0.7, merge_dist=None)




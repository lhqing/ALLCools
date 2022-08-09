:py:mod:`ALLCools.reptile.reptile`
==================================

.. py:module:: ALLCools.reptile.reptile


Module Contents
---------------

.. py:function:: _create_train_region_ds(reptile)


.. py:function:: _create_train_dmr_ds(reptile, train_regions_bed, train_label)


.. py:function:: _create_query_region_ds(reptile)


.. py:function:: _create_query_dmr_ds(reptile, dmr_regions_bed_df)


.. py:function:: _get_data_and_label(region_ds, modalities, sample, fillna_by_zero_list)


.. py:function:: _predict_sample(region_ds_path, region_dim, modalities, fillna_by_zero, sample, output_path, mask_cutoff=0.3, chunk_size=100000)


.. py:function:: _call_enhancer_region(bw_path, dmr_bed, threshold, merge_dist, chrom_size_path)


.. py:class:: REPTILE(output_path, train_regions, dmr_regions, train_region_labels, train_sample, bigwig_table, chrom_size_path, window_size=2000, step_size=200, dmr_slop=150, fillna_by_zero=None)

   REPTILE Pipeline for Enhancer Prediction

   .. py:method:: generate_region_ds()

      Generate RegionDS for training and query.


   .. py:method:: train_region_ds()
      :property:

      Get train region RegionDS.


   .. py:method:: train_dmr_ds()
      :property:

      Get train DMR RegionDS.


   .. py:method:: query_region_ds()
      :property:

      Get query region RegionDS.


   .. py:method:: query_dmr_ds()
      :property:

      Get query DMR RegionDS.


   .. py:method:: region_model()
      :property:

      Get the region model.


   .. py:method:: dmr_model()
      :property:

      Get the DMR model.


   .. py:method:: _validate_region_name(name)


   .. py:method:: annotate_by_bigwigs(name, slop, cpu, redo=False)

      Annotate genome regions by bigwigs.


   .. py:method:: _filter_na_train(name, sample, max_na_rate=0.5)


   .. py:method:: prepare_training_input(name)

      Prepare training input for a type of region.


   .. py:method:: auto_ml(data, label, output_path, train_size=0.75, random_state=42, cpu=1, tpot_generations=5, tpot_max_time_mins=60, **tpot_kwargs)
      :staticmethod:

      Perform the auto-ML training and save the model to output_path.


   .. py:method:: _train(region_dim, slop, cpu, **kwargs)


   .. py:method:: train_region_model(slop=None, cpu=1, **kwargs)

      Train a model for genomic regions.


   .. py:method:: train_dmr_model(slop=None, cpu=1, **kwargs)

      Train a model for DMRs.


   .. py:method:: fit(cpu=10, **kwargs)

      Train everything by default parameters.


   .. py:method:: _predict(region_dim, cpu, mask_cutoff)


   .. py:method:: predict(cpu, mask_cutoff=0.3, bw_bin_size=10, enhancer_cutoff=0.7)

      Predict enhancer from score tracks.


   .. py:method:: _dump_sample(sample, mask_cutoff, bw_bin_size)


   .. py:method:: dump_bigwigs(cpu, mask_cutoff, bw_bin_size)

      Dump bigwig files for each sample.


   .. py:method:: call_enhancers(threshold=0.7, merge_dist=None)

      Call enhancers from REPTILE score tracks.




:py:mod:`ALLCools.sandbox.motif.utilities`
==========================================

.. py:module:: ALLCools.sandbox.motif.utilities


Module Contents
---------------

.. py:function:: meme_to_homer(meme_path, homer_path, score_power=0.85)

   Transfer MEME motif format into Homer motif format.
   Based on description here: http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html
   The score_power controls Log odds detection threshold by max_score ** score_power

   :param meme_path: Input path in meme format
   :param homer_path: Output path in homer format
   :param score_power: Log odds detection threshold = max_score ** score_power


.. py:function:: meme_motif_file_to_dict(meme_motif_paths)


.. py:function:: single_meme_txt_to_pfm_df(text, bits_scale=True)


.. py:function:: meme_to_pfm_dict(meme_motif_paths, bits_scale=True)


.. py:function:: plot_pfm(pfm, ax=None, logo_kws=None)


.. py:function:: split_meme_motif_file(meme_motif_paths, output_dir)

   Given multi motif meme format file, split into single motif meme format file

   :param meme_motif_paths:
   :param output_dir:



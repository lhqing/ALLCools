:py:mod:`ALLCools.count_matrix.snap`
====================================

.. py:module:: ALLCools.count_matrix.snap


Module Contents
---------------

.. py:data:: SNAP_BD_FIELDS
   

   

.. py:function:: read_snap(file_path, bin_kind: Union[int, str] = 5000)

   Read snap hdf5 file into anndata.AnnData.


.. py:function:: _read_snap_meta(f)


.. py:function:: _read_snap_genes(file_path)

   Read gene data from snap hdf5 file into anndata.AnnData.


.. py:function:: _read_snap_bins(file_path, bin_size=5000)

   Read bin data from snap hdf5 file into anndata.AnnData


.. py:function:: adata_to_df(adata, var_dim, obs_dim='cell', dtype=np.uint8)


.. py:function:: snap_to_xarray(snap_path, bin_sizes=(5000, ), gene=True, dtype=np.uint8)


.. py:function:: snap_to_zarr(snap_path, output_path, bin_sizes=(5000, ), gene=True, dtype=np.uint8, index_prefix=None)



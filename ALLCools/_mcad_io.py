"""
The hdf5 functions copied from AnnData
https://github.com/theislab/anndata

BSD 3-Clause License

Copyright (c) 2017-2018 P. Angerer, F. Alexander Wolf, Theis Lab
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import pathlib
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as ss
from anndata import AnnData
from anndata import h5py
from collections.abc import Mapping
from pandas.api.types import is_string_dtype, is_categorical
from scipy.sparse import issparse

from .utilities import generate_chrom_bin_bed_dataframe


def df_to_records_fixed_width(df):
    """
    Turn pandas DataFrame into numpy recarray, and take care of data types
    """
    uns = {}  # unstructured dictionary for storing categories
    names = ['index']
    if is_string_dtype(df.index):
        max_len_index = df.index.map(len).max()
        index = df.index.values.astype('S{}'.format(max_len_index))
    else:
        index = df.index.values
    arrays = [index]
    for k in df.columns:
        names.append(k)
        if is_string_dtype(df[k]):
            c = df[k].astype('category')
            # transform to category only if there are less categories than
            # observations
            if len(c.cat.categories) < len(c):
                uns[k + '_categories'] = c.cat.categories.values
                arrays.append(c.cat.codes)
            else:
                max_len_index = df[k].map(len).max()
                arrays.append(df[k].values.astype('S{}'.format(max_len_index)))
        elif is_categorical(df[k]):
            uns[k + '_categories'] = df[k].cat.categories
            arrays.append(df[k].cat.codes)
        else:
            arrays.append(df[k].values)
    formats = [v.dtype for v in arrays]
    return np.rec.fromarrays(
        arrays,
        dtype={'names': names, 'formats': formats}), uns


def _to_dict_fixed_width_arrays(adata):
    """A dict of arrays that stores data and annotation.

    It is sufficient for reconstructing the object.
    """
    obs_rec, uns_obs = df_to_records_fixed_width(adata.obs)
    var_rec, uns_var = df_to_records_fixed_width(adata.var)
    d = {
        'X': adata.X,
        'obs': obs_rec,
        'var': var_rec,
        'obsm': adata.obsm,
        'varm': adata.varm,
        'layers': adata,
        # add the categories to the unstructured annotation
        'uns': {**adata.uns, **uns_obs, **uns_var}}
    return d


def _write_mcad(file_path,
                key_base,
                adata: AnnData,
                force_dense: bool = False,
                **kwargs):
    file_path = pathlib.Path(file_path)
    if not file_path.name.endswith('.mcad'):
        raise ValueError("Filename needs to end with '.mcad'.")
    if adata.isbacked:
        # close so that we can reopen below
        adata.file.close()
    d = _to_dict_fixed_width_arrays(adata)
    # we're writing to a different location than the backing file
    # - load the matrix into the memory...
    if adata.isbacked and file_path != adata.filename:
        d['X'] = adata.X[:]
    # need to use 'a' if backed, otherwise we loose the backed objects
    with h5py.File(file_path, 'a', force_dense=force_dense) as f:
        for key, value in d.items():
            key = key_base + '/' + key
            _write_key_value_to_h5(f, key, value, **kwargs)


def _preprocess_writing(value):
    if value is None or issparse(value):
        return value
    else:
        value = np.array(value)  # make sure value is an array
        if value.ndim == 0:
            value = np.array([value])  # hm, why that?
    # make sure string format is chosen correctly
    if value.dtype.kind == 'U':
        value = value.astype(h5py.special_dtype(vlen=str))
    return value


def _write_key_value_to_h5(f, key, value, **kwargs):
    if isinstance(value, Mapping):
        for k, v in value.items():
            if not isinstance(k, str):
                warnings.warn(
                    'dict key {} transformed to str upon writing to h5,'
                    'using string keys is recommended'
                        .format(k)
                )
            _write_key_value_to_h5(f, key + '/' + str(k), v, **kwargs)
        return

    value = _preprocess_writing(value)

    # for some reason, we need the following for writing string arrays
    if key in f.keys() and value is not None:
        del f[key]

    # ignore arrays with empty dtypes
    if value is None or not value.dtype.descr:
        return
    try:
        if key in set(f.keys()):
            is_valid_group = (
                    isinstance(f[key], h5py.Group)
                    and f[key].shape == value.shape
                    and f[key].dtype == value.dtype
                    and not isinstance(f[key], h5py.SparseDataset)
            )
            if not is_valid_group and not issparse(value):
                f[key][()] = value
                return
            else:
                del f[key]
        f.create_dataset(key, data=value, **kwargs)
    except TypeError:
        try:
            if value.dtype.names is None:
                dt = h5py.special_dtype(vlen=str)
                f.create_dataset(key, data=value, dtype=dt, **kwargs)
            else:  # try writing composite datatypes with byte strings
                new_dtype = [
                    (dt[0], 'S{}'.format(int(dt[1][2:]) * 4))
                    for dt in value.dtype.descr
                ]
                if key in set(f.keys()):
                    if (
                            f[key].shape == value.shape
                            and f[key].dtype == value.dtype
                    ):
                        f[key][()] = value.astype(new_dtype)
                        return
                    else:
                        del f[key]
                f.create_dataset(
                    key, data=value.astype(new_dtype), **kwargs)
        except Exception as e:
            warnings.warn(
                'Could not save field with key = {!r} to hdf5 file: {}'
                    .format(key, e)
            )


def csr_matrix_to_mcad(key_base, matrix_paths, obs_names, chrom_size_path, bin_size, output_path,
                       compression=None, compression_opts=None):
    """
    Save a list of csc matrix into specific place (key_base) of the hdf5 file (output_path)
    Not doing param check here, only call this function from aggregate_region_count_to_mcad
    """

    if pathlib.Path(output_path).exists():
        return output_path

    var_df = generate_chrom_bin_bed_dataframe(chrom_size_path=chrom_size_path,
                                              window_size=bin_size)
    total_matrix = ss.vstack([ss.load_npz(path) for path in matrix_paths])

    adata = AnnData(X=total_matrix,
                    obs=pd.DataFrame([], index=obs_names),
                    var=var_df[['chrom']],
                    uns=dict(bin_size=bin_size,
                             chrom_size_path=chrom_size_path))
    _write_mcad(file_path=output_path,
                key_base=key_base,
                adata=adata,
                force_dense=False,
                compression=compression,
                compression_opts=compression_opts)
    return output_path

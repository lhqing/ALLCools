import pathlib
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import permutations
from typing import Tuple

import numpy as np
import pandas as pd
import pysam
import zarr
from pysam import TabixFile
from scipy.sparse import csc_matrix

from ALLCools.utilities import reverse_complement


def _get_dtype_max(dtype):
    try:
        max_value = np.iinfo(dtype).max
    except ValueError:
        max_value = np.finfo(dtype).max
    return max_value


class _ALLCFile(TabixFile):
    def __init__(self, filename):
        # validate path and tbi file
        self._path = pathlib.Path(filename)
        self._path_tbi = pathlib.Path(f"{self._path}.tbi")
        if not self._path.exists():
            raise FileNotFoundError(f"{self._path} not found")
        if not self._path_tbi.exists():
            raise FileNotFoundError(f"Index file {self._path_tbi} not found")

        # noinspection PyArgumentList
        super().__init__()
        self.handle = None
        return

    def fetch(self, *args, **kwargs):
        if self.closed:
            raise RuntimeError("File not open")

        row_iter = super().fetch(*args, **kwargs)
        for row in row_iter:
            # chrom, pos, strand, context, mc, cov, p
            ll = row.strip().split("\t")
            ll[1] = int(ll[1])
            ll[4] = int(ll[4])
            ll[5] = int(ll[5])
            yield ll

    def fetch_df(self, *args, **kwargs) -> pd.DataFrame:
        """
        Fetch ALLC files by position and return a pandas dataframe.

        Parameters
        ----------
        args :
            positional arguments for `TabixFile.fetch`
        kwargs :
            keyword arguments for `TabixFile.fetch`

        Returns
        -------
        df :
            a pandas dataframe with shape (n_rows, 7).
        """
        df = pd.DataFrame(
            list(self.fetch(*args, **kwargs)),
            columns=["chrom", "pos", "strand", "context", "mc", "cov", "p"],
        )
        return df

    def fetch_pos_and_count_array(self, *args, pos_base=0, **kwargs) -> np.ndarray:
        """
        Fetch ALLC files by position and return a numpy array with shape (n_rows, 3).

        The first column is the position.
        If pos_base is 0, the position is converted to 0-base.Original ALLC files are 1-base.
        the second column is the mc count.
        the third column is the cov count.

        Parameters
        ----------
        args :
            positional arguments for `TabixFile.fetch`
        pos_base :
            the position base for the array, default is 0, which means the array pos column is 0-based.
            ALLC file is 1-based.
        kwargs :
            keyword arguments for TabixFile.fetch

        Returns
        -------
        records :
            a numpy array with shape (n_rows, 3).
        """
        try:
            records = np.array([(pos, mc, cov) for _, pos, _, _, mc, cov, _ in self.fetch(*args, **kwargs)])
        except ValueError as e:
            print(f"{self.filename} fetch causing ValueError")
            raise e
        if records.shape[0] == 0:
            # no records in allc
            # this records is still empty, but shape is correct
            records = np.zeros(shape=(0, 3))

        if pos_base == 0:
            # ALLC pos is 1-based,
            # the pos 1 actually point to base 0 in 0-based index.
            records[:, 0] -= 1

        return records

    def fetch_full_base_array(
        self, *args, pos_min=None, pos_max=None, sparse=True, dtype="float32", **kwargs
    ) -> Tuple[np.ndarray, int]:
        """
        Fetch ALLC files by position and return a full-base numpy array with shape (pos_max, 2).

        The first column is the mc count. The second column is cov count.

        Parameters
        ----------
        args :
            positional arguments for `TabixFile.fetch`
        pos_min :
            the minimum position fetched from ALLC or the minimum position in the query region provided by user.
            Use this number as the offset for the array position.
        pos_max :
            the relative maximum position (compare to pos_min) fetched from ALLC
            or the maximum position in the query region provided by user.
        sparse :
            whether to return a sparse array or a dense array.
        dtype :
            the data type of the array.
        kwargs :
            keyword arguments for TabixFile.fetch

        Returns
        -------
        full_base_array :
            a full-base numpy array with shape (pos_max, 2).
            If no data fetched from ALLC, return empty array with shape (0, 2)
        """
        # pos_base has to be 0
        kwargs["pos_base"] = 0
        _a = self.fetch_pos_and_count_array(*args, **kwargs)
        if _a.shape[0] == 0:
            return np.zeros(shape=(0, 2)), pos_min

        pos = _a[:, 0]
        mc = _a[:, 1]
        cov = _a[:, 2]

        # prevent dtype overflow at extremely high coverage sites
        dtype_max = _get_dtype_max(dtype)
        cov_max = cov.max()
        if cov_max > dtype_max:
            # warning about potential data overflow
            import warnings

            warnings.warn(
                f"cov max {cov_max} is larger than dtype max {dtype_max}, "
                "in order to prevent overflow, mc and cov sites with cov > dtype_max will be scale down."
            )

            # scale down mc and cov value if cov is larger than dtype_max
            _flag = cov > dtype_max
            cov[_flag] = np.round(cov[_flag] / cov_max * dtype_max).astype(cov.dtype)
            mc[_flag] = np.round(mc[_flag] / cov_max * dtype_max).astype(mc.dtype)

        # pos_min use to determine the offset of the array position
        _pos_min = pos.min()
        if pos_min is None:
            pos_min = _pos_min
        else:
            if pos_min > _pos_min:
                raise ValueError(f"pos_min {pos_min} is larger than the minimum position {_pos_min}")

        # apply offset to the array position
        pos -= pos_min

        # pos_max use to determine shape of full_base_array
        _pos_max = pos.max()
        if pos_max is None:
            pos_max = _pos_max
        else:
            if pos_max < _pos_max:
                raise ValueError(f"pos_max {pos_max} is smaller than the maximum position {_pos_max}")

        # save the mc and cov count in a full base sparse array
        full_base_array = csc_matrix(
            (np.concatenate([mc, cov]), np.concatenate([pos, pos]), np.cumsum([0, pos.size, pos.size])),
            shape=(pos_max, 2),
            dtype=dtype,
        )

        if not sparse:
            # turn to numpy array
            full_base_array = full_base_array.toarray()
        return full_base_array, pos_min


class _BaseDSChunkWriter:
    """Base class for writing chunks to a dataset."""

    def __init__(self, zarr_path, sample_chunk: pd.Series, sample_offset: int, chrom_chunk: list, dtype):
        self.zarr_path = zarr_path
        self.sample_chunk = sample_chunk
        self.sample_offset = sample_offset
        self.n_sample = sample_chunk.size

        self.chrom, self.start, self.end = chrom_chunk

        self.dtype = dtype

        self.execute()
        return

    def _get_data(self):
        pos_min = self.start
        pos_max = self.end - self.start

        # shape (pos, sample, count_type)
        _data = np.zeros(shape=(pos_max, self.n_sample, 2), dtype=self.dtype)
        for i, path in enumerate(self.sample_chunk.values):
            with _ALLCFile(path) as f:
                # read region count as a full_base_array from each allc file and save in _data
                full_base_array, _ = f.fetch_full_base_array(
                    reference=self.chrom,
                    start=self.start,
                    end=self.end,
                    pos_min=pos_min,
                    pos_max=pos_max,
                    sparse=False,
                    dtype=self.dtype,
                )
                if full_base_array.shape[0] == 0:
                    # no data
                    continue
                _data[:, i, :] += full_base_array
        return _data

    def execute(self):
        # zarr_da shape (pos, sample, count_type)
        zarr_da = zarr.open(self.zarr_path, mode="r+")

        pos_slice = slice(self.start, self.end)
        sample_slice = slice(self.sample_offset, self.sample_offset + self.n_sample)
        zarr_da[pos_slice, sample_slice, :] = self._get_data()
        return


def _create_contexts(context_size, c_pos, include_n):
    if include_n:
        bases = "ATCGN"
    else:
        bases = "ATCG"

    contexts = sorted({"".join(c) for c in permutations(bases * context_size, context_size) if c[c_pos] == "C"})
    return bases, contexts


def _create_codebook_single_chrom(zarr_root, chrom, genome_fasta_path, chrom_size, chunk_size, contexts):
    print(f"Generating codebook for {chrom}...")
    with pysam.FastaFile(genome_fasta_path, filepath_index=f"{genome_fasta_path}.fai") as fa:
        sequence = fa.fetch(chrom).upper()

    # make sure chrom size consistent
    assert len(sequence) == chrom_size, f"chrom {chrom} size {len(sequence)} is not consistent with {chrom_size}"

    context_codebook = zarr_root.require_dataset(
        "codebook", shape=(chrom_size, len(contexts)), chunks=chunk_size, dtype="int8"
    )
    context_codebook.attrs["_ARRAY_DIMENSIONS"] = ["pos", "mc_type"]

    for i, context in enumerate(contexts):
        print(context, end=", ")
        # positive strand
        # in order to count overlapping contexts, this regex pattern is necessary
        # see https://stackoverflow.com/questions/5616822/how-to-use-regex-to-find-all-overlapping-matches
        pos_context_str = f"(?=({context}))"
        pos_pattern = re.compile(pos_context_str)
        pos_pos = [match.start() for match in pos_pattern.finditer(sequence)]
        # negative strand
        neg_context_str = f"(?=({reverse_complement(context)}))"
        neg_pattern = re.compile(neg_context_str)
        neg_dodge = len(context) - 1
        neg_pos = [match.end() + neg_dodge for match in neg_pattern.finditer(sequence)]

        pos_value = pd.Series({i: 1 for i in pos_pos})
        neg_value = pd.Series({i: -1 for i in neg_pos})
        value = pd.concat([pos_value, neg_value]).sort_index()
        context_codebook[value.index, i] = value.values

    context_codebook.attrs["contexts"] = list(contexts)
    print()
    return


def _create_mc_type_zarr(root, contexts, c_pos, include_n, context_size):
    n_type = len(contexts)
    mc_type = root.require_dataset("mc_type", shape=(n_type,), chunks=(n_type,), dtype="<U50")
    mc_type[:] = np.array(contexts)
    mc_type.attrs["_ARRAY_DIMENSIONS"] = "mc_type"
    mc_type.attrs["c_pos"] = c_pos
    mc_type.attrs["include_n"] = include_n
    mc_type.attrs["context_size"] = context_size
    return


class _CreateSingleChromosomeBaseDS:
    """Create the BaseDS zarr dataset for one chromosome."""

    def __init__(
        self,
        path,
        chrom,
        chrom_size,
        sample_path_series: pd.Series,
        fasta_path=None,
        generate_codebook=True,
        context_size=3,
        c_pos=0,
        include_n=False,
        pos_chunk=10000000,
        sample_chunk=10,
        data_dtype="i4",
        mode="a",
        n_cpu=1,
    ):
        self.path = path
        self.data_path = f"{self.path}/data"
        self.root = zarr.open(path, mode=mode)
        self.pos_chunk = pos_chunk
        self.sample_chunk = sample_chunk
        self.n_cpu = n_cpu
        self.data_dtype = data_dtype

        self.log_dir_path = pathlib.Path(path) / ".log"
        self.log_dir_path.mkdir(exist_ok=True)

        # chrom size info
        self.chrom = chrom
        self.chrom_size = chrom_size

        # mC type context info
        self.c_pos = c_pos
        self.context_size = context_size
        self.include_n = include_n
        self.bases, self.contexts = self._get_contexts()

        # chrom chunk info
        self._generate_codebook = generate_codebook
        self._genome_fasta_path = fasta_path
        if self._generate_codebook:
            assert self._genome_fasta_path is not None, "genome fasta path is required when generating codebook"
        self.chrom_chunk_pos = self._get_chrom_chunk_pos()

        # sample info
        self.sample_path_series = sample_path_series
        self.sample_ids = self.sample_path_series.index.tolist()
        self.n_sample = len(self.sample_ids)

        # get chrom pos and sample chunks for iterating writing
        self.chrom_chunks = self._get_chrom_chunks()
        self.sample_chunks = self._get_sample_chunks()

        # init zarr structure
        self._init_zarr()

        # write data in parallel
        self._write_data()
        return

    def _get_contexts(self):
        bases, contexts = _create_contexts(context_size=self.context_size, c_pos=self.c_pos, include_n=self.include_n)
        return bases, contexts

    def _create_contexts(self):
        bases, contexts = self._get_contexts()
        _create_mc_type_zarr(
            root=self.root,
            contexts=contexts,
            c_pos=self.c_pos,
            include_n=self.include_n,
            context_size=self.context_size,
        )
        return bases, contexts

    def _create_codebook(self):
        _create_codebook_single_chrom(
            zarr_root=self.root,
            chrom=self.chrom,
            genome_fasta_path=self._genome_fasta_path,
            chrom_size=self.chrom_size,
            chunk_size=(self.pos_chunk, len(self.bases)),
            contexts=self.contexts,
        )
        return

    def _get_chrom_chunk_pos(self):
        # create chunk pos array
        chrom_chunk_pos = list(range(0, self.chrom_size, self.pos_chunk))
        if chrom_chunk_pos[-1] != self.chrom_size:
            chrom_chunk_pos.append(self.chrom_size)
        chunk_pos = self.root.require_dataset(
            "chunk_pos", shape=(len(chrom_chunk_pos),), chunks=(len(chrom_chunk_pos),), dtype="<i4"
        )
        chunk_pos[:] = chrom_chunk_pos
        chunk_pos.attrs["_ARRAY_DIMENSIONS"] = "chunk_pos"
        return chrom_chunk_pos

    def _get_chrom_chunks(self):
        chrom_chunks = []
        for i in range(len(self.chrom_chunk_pos) - 1):
            start = self.chrom_chunk_pos[i]
            end = self.chrom_chunk_pos[i + 1]
            chrom_chunks.append([self.chrom, start, end])
        return chrom_chunks

    def _get_sample_chunks(self):
        # create chunk pos array
        sample_chunks = []
        for chunk_start in range(0, self.n_sample, self.sample_chunk):
            sample_chunk = self.sample_path_series[chunk_start : chunk_start + self.sample_chunk]
            sample_chunks.append(sample_chunk)
        return sample_chunks

    def _create_count_type(self):
        count_type = self.root.require_dataset("count_type", shape=(2,), chunks=(2,), dtype="<U10")
        count_type[:] = ["mc", "cov"]
        count_type.attrs["_ARRAY_DIMENSIONS"] = ["count_type"]
        return

    def _create_sample_id(self):
        sample_id = self.root.require_dataset(
            "sample_id", shape=(self.n_sample,), chunks=(self.n_sample,), dtype="<U256"
        )
        sample_id[:] = self.sample_ids
        sample_id.attrs["_ARRAY_DIMENSIONS"] = "sample_id"
        return

    def _create_data(self):
        data_name = pathlib.Path(self.data_path).name
        z = self.root.require_dataset(
            data_name,
            shape=(self.chrom_size, self.n_sample, 2),
            chunks=(self.pos_chunk, self.sample_chunk, 1),
            dtype=self.data_dtype,
        )
        z.attrs["_ARRAY_DIMENSIONS"] = ["pos", "sample_id", "count_type"]
        return

    def _add_root_attrs(self):
        self.root.attrs["context_size"] = self.context_size
        self.root.attrs["c_pos"] = self.c_pos
        self.root.attrs["include_n"] = self.include_n

        self.root.attrs["data_dtype"] = self.data_dtype
        self.root.attrs["obs_dim"] = "sample_id"
        self.root.attrs["obs_size"] = self.n_sample
        self.root.attrs["obs_chunk"] = self.sample_chunk

        self.root.attrs["chrom"] = self.chrom
        self.root.attrs["chrom_size"] = self.chrom_size
        self.root.attrs["pos_chunk"] = self.pos_chunk
        return

    def _consolidate_metadata(self):
        zarr.consolidate_metadata(self.path)
        return

    def _init_zarr(self):
        success_flag = self.log_dir_path / "_PREPARE_SUCCESS"
        if success_flag.exists():
            print("Zarr structure already initiated. Skipping...")
            return

        # create the zarr dataset structure
        self._create_contexts()
        if self._generate_codebook:
            self._create_codebook()
        self._create_count_type()
        self._create_sample_id()
        # init empty sample data array per chromosome
        self._create_data()

        # add attrs
        self._add_root_attrs()

        # consolidate metadata
        self._consolidate_metadata()

        # touch success flag
        success_flag.touch()
        return

    def _write_data(self):
        """Write data to zarr in parallel."""
        # config blosc compressor when run multi-processing
        # see zarr doc here: https://zarr.readthedocs.io/en/stable/tutorial.html#configuring-blosc
        if self.n_cpu > 1:
            from numcodecs import blosc

            blosc.use_threads = False

        final_success_flag = self.log_dir_path / "_WRITE_SUCCESS"
        if final_success_flag.exists():
            print("Data already written. Skipping...")
            return

        with ProcessPoolExecutor(max_workers=self.n_cpu) as executor:
            futures = {}
            for i, chrom_chunk in enumerate(self.chrom_chunks):
                sample_offset = 0
                for j, sample_chunk in enumerate(self.sample_chunks):
                    success_flag = self.log_dir_path / f"_WRITE_SUCCESS_{i}_{j}"
                    if success_flag.exists():
                        print(
                            f"BaseDS Writer: {self.chrom} region {i}/{len(self.chrom_chunks)} "
                            f"sample {j}/{len(self.sample_chunks)} already exist."
                        )
                        continue
                    f = executor.submit(
                        _BaseDSChunkWriter,
                        zarr_path=str(self.data_path),
                        chrom_chunk=chrom_chunk,
                        sample_offset=sample_offset,
                        sample_chunk=sample_chunk,
                        dtype=self.data_dtype,
                    )
                    sample_offset += len(sample_chunk)
                    futures[f] = (i, j)

            for future in as_completed(futures):
                i, j = futures[future]
                print(
                    f"BaseDS Writer: {self.chrom} region {i}/{len(self.chrom_chunks)} "
                    f"sample {j}/{len(self.sample_chunks)} done."
                )
                future.result()

                # create a tag for successful completion
                success_flag = self.log_dir_path / f"_WRITE_SUCCESS_{i}_{j}"
                success_flag.touch()

        # create final success flag
        final_success_flag.touch()

        # cleanup individual success flags
        for success_flag in self.log_dir_path.glob("_WRITE_SUCCESS_*"):
            success_flag.unlink()
        return


def generate_base_ds(
    output_dir_path,
    genome_fasta_path,
    chrom_size_path,
    allc_table_path,
    generate_codebook=True,
    context_size=3,
    c_pos=0,
    include_n=False,
    pos_chunk=1000000,
    sample_chunk=200,
    data_dtype="uint16",
    mode="a",
    n_cpu=1,
):
    """Create a BaseDS zarr dataset."""
    output_dir_path = pathlib.Path(output_dir_path)
    output_dir_path.mkdir(exist_ok=True, parents=True)

    chroms = pd.read_csv(chrom_size_path, header=None, sep="\t", index_col=0, names=["chrom", "size"])["size"]
    # save chroms to output_dir_path
    chroms.to_csv(f"{output_dir_path}/chrom_sizes.csv", index=True, header=False)

    allc_table_path = pathlib.Path(allc_table_path)
    sample_path_series = pd.read_csv(
        allc_table_path,
        header=None,
        index_col=0,
        sep="\t" if allc_table_path.name.endswith("tsv") else ",",
    ).squeeze()

    # check if all the allc paths in sample_path_series exist
    flag = False
    for p in sample_path_series.to_list():
        if not pathlib.Path(p).exists():
            flag = True
            print(f"{p} does not exist.")
    if flag:
        raise FileNotFoundError(f"Some of the allc files in {allc_table_path} do not exist.")

    if generate_codebook:
        fasta_path = pathlib.Path(genome_fasta_path)
        fai_path = pathlib.Path(f"{fasta_path}.fai")
        if not fai_path.exists():
            print(f"{fai_path} not found. Try generating...")
            import subprocess

            try:
                subprocess.run(["samtools", "faidx", str(fasta_path)], capture_output=True, check=True)
            except subprocess.CalledProcessError as e:
                print(
                    f"samtools faidx failed with error: {e.stderr.decode()}. "
                    f"Please try to generate the fai file manually."
                )
                raise e
    else:
        fasta_path = None

    output_dir_path = pathlib.Path(output_dir_path)
    output_dir_path.mkdir(exist_ok=True, parents=True)

    for chrom, chrom_size in chroms.items():
        chrom_dir_path = f"{output_dir_path}/{chrom}"
        print(f"Creating {chrom}, size: {chrom_size}")
        _CreateSingleChromosomeBaseDS(
            path=chrom_dir_path,
            chrom=chrom,
            chrom_size=chrom_size,
            fasta_path=fasta_path,
            sample_path_series=sample_path_series,
            generate_codebook=generate_codebook,
            context_size=context_size,
            c_pos=c_pos,
            include_n=include_n,
            pos_chunk=pos_chunk,
            sample_chunk=sample_chunk,
            data_dtype=data_dtype,
            mode=mode,
            n_cpu=n_cpu,
        )
    return


def generate_mc_codebook(
    path,
    genome_fasta_path,
    chrom_size_path,
    mode="w",
    context_size=3,
    c_pos=0,
    include_n=False,
    pos_chunk=1000000,
    n_cpu=1,
):
    """Generate mC context codebook for BaseDS."""
    bases, contexts = _create_contexts(context_size=context_size, c_pos=c_pos, include_n=include_n)
    chrom_sizes = pd.read_csv(chrom_size_path, header=None, sep="\t", index_col=0, names=["chrom", "size"]).squeeze()

    if n_cpu > 1:
        from numcodecs import blosc

        blosc.use_threads = False

    with ProcessPoolExecutor(max_workers=n_cpu) as executor:
        futures = {}
        for chrom, chrom_size in chrom_sizes.items():
            chrom_path = pathlib.Path(f"{path}/{chrom}")
            chrom_path.mkdir(exist_ok=True, parents=True)
            chrom_root = zarr.open(str(chrom_path), mode=mode)

            future = executor.submit(
                _create_codebook_single_chrom,
                zarr_root=chrom_root,
                chrom=chrom,
                genome_fasta_path=genome_fasta_path,
                chrom_size=chrom_size,
                chunk_size=(pos_chunk, len(bases)),
                contexts=contexts,
            )
            futures[future] = chrom_root, chrom_path

        for future in as_completed(futures):
            chrom_root, chrom_path = futures[future]
            # add mc type index
            _create_mc_type_zarr(
                root=chrom_root, contexts=contexts, c_pos=c_pos, include_n=include_n, context_size=context_size
            )
            zarr.consolidate_metadata(str(chrom_path))
    return

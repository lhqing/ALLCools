import re
from itertools import permutations

import numpy as np
import zarr


class _BaseDSChunkWriter:
    """Base class for writing chunks to a dataset."""

    def __init__(self):
        pass


class _CreateSingleChromosomeBaseDS:
    """Create the BaseDS zarr dataset for one chromosome."""

    def __init__(
        self,
        path,
        chrom,
        chrom_size,
        chrom_sequence,
        sample_ids,
        context_size=3,
        c_pos=0,
        include_n=False,
        pos_chunk=10000000,
        sample_chunk=10,
        data_dtype="i4",
    ):
        self.path = path
        self.root = zarr.open(path, mode="a")
        self.pos_chunk = pos_chunk
        self.sample_chunk = sample_chunk

        # chrom size info
        self.chrom = chrom
        self.chrom_size = chrom_size

        # mC type context info
        self.c_pos = c_pos
        self.context_size = context_size
        self.include_n = include_n
        self.bases, self.contexts = self._create_contexts()

        # codebook per chromosome
        self.chrom_sequence = chrom_sequence.upper()
        self._create_codebook()

        # count type
        self._create_count_type()

        # sample info
        self.sample_ids = sample_ids
        self.n_sample = len(sample_ids)
        self._create_sample_id()

        # init empty sample dataarray per chromosome
        self.data_dtype = data_dtype
        self._create_data()

        # consolidate metadata
        self._consolidate_metadata()

        # write data in parallel
        self._write_data()
        return

    def _get_contexts(self):
        if self.include_n:
            bases = "ATCGN"
        else:
            bases = "ATCG"

        contexts = sorted(
            {"".join(c) for c in permutations(bases * self.context_size, self.context_size) if c[self.c_pos] == "C"}
        )
        return bases, contexts

    def _create_contexts(self):
        bases, contexts = self._get_contexts()
        n_type = len(contexts)
        mc_type = self.root.require_dataset("mc_type", shape=(n_type,), chunks=(n_type,), dtype="<U50")
        mc_type[:] = np.array(contexts)
        mc_type.attrs["_ARRAY_DIMENSIONS"] = "mc_type"
        mc_type.attrs["c_pos"] = self.c_pos
        mc_type.attrs["include_n"] = self.include_n
        mc_type.attrs["context_size"] = self.context_size
        return bases, contexts

    def _create_codebook(self):
        chrom_sequence = self.chrom_sequence

        contexts = self.contexts
        print(f"Generating codebook for {self.chrom}...")

        # make sure chrom size consistent
        chrom_size = self.chrom_size
        assert len(chrom_sequence) == chrom_size

        context_codebook = self.root.require_dataset(
            "codebook", shape=(chrom_size, len(contexts)), chunks=(self.pos_chunk, len(self.bases)), dtype="bool"
        )
        context_codebook.attrs["_ARRAY_DIMENSIONS"] = ["pos", "mc_type"]

        for i, context in enumerate(contexts):
            pattern = re.compile(context)
            context_pos = [match.start() for match in pattern.finditer(chrom_sequence)]
            context_codebook[context_pos, i] = True

        # create chunk pos array
        chrom_chunk_pos = list(range(0, chrom_size, self.pos_chunk))
        if chrom_chunk_pos[-1] != chrom_size:
            chrom_chunk_pos.append(chrom_size)
        chunk_pos = self.root.require_dataset(
            "chunk_pos", shape=(len(chrom_chunk_pos),), chunks=(len(chrom_chunk_pos),), dtype="<i4"
        )
        chunk_pos[:] = chrom_chunk_pos
        chunk_pos.attrs["_ARRAY_DIMENSIONS"] = "chunk_pos"
        return

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
        z = self.root.require_dataset(
            "data",
            shape=(self.chrom_size, self.n_sample, len(self.contexts)),
            chunks=(self.pos_chunk, self.sample_chunk, 1),
            dtype=self.data_dtype,
        )
        z.attrs["_ARRAY_DIMENSIONS"] = ["pos", "sample_id", "count_type"]
        return

    def _consolidate_metadata(self):
        zarr.consolidate_metadata(self.path)
        return

    def _write_data(self):
        """Parallel write data to zarr."""
        # determine chunks
        # process each chunk in parallel
        # using BaseDSChunkWriter class
        return


class CreateSingleChromosomeBaseDS:
    def __init__(self):
        pass

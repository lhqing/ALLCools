from Bio import motifs, SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np
import xarray as xr
import pathlib
import subprocess
from ..mcds.utilities import write_ordered_chunks
from ..mcds import RegionDS


def _get_motif_threshold(motif, method, threshold_value):
    if method == "patser":
        value = motif.pssm.distribution().threshold_patser()
    elif method == "balanced":
        value = motif.pssm.distribution().threshold_balanced(threshold_value)
    elif method == "fpr":
        value = motif.pssm.distribution().threshold_fpr(threshold_value)
    elif method == "fnr":
        value = motif.pssm.distribution().threshold_fnr(threshold_value)
    else:
        raise ValueError(
            f'Method needs to be ["patser", "balanced", "fpr", "fnr"], got {method}. '
            f"Check this page for more details: https://biopython-tutorial.readthedocs.io/en/latest/"
            f"notebooks/14%20-%20Sequence%20motif%20analysis%20using%20Bio.motifs.html"
            f"#Selecting-a-score-threshold"
        )
    return value


def _single_seq_single_motif_max_score(seq, motif, threshold):
    if len(seq) < motif.length:
        return 0, 0, 0
    try:
        # position and direction information can be collected in motif.pssm.search
        # but not used here
        scores = [i[1] for i in motif.pssm.search(seq, threshold=threshold)]
        n_sites = len(scores)
        max_score = max(scores)
        total_score = sum(scores)
        return n_sites, max_score, total_score
    except ValueError:
        # motif not found
        return 0, 0, 0


class MotifSet:
    def __init__(self, motifs, meta_table=None, motif_cluster_col="cluster_id"):
        self.motif_list = motifs
        self.motif_dict = {motif.name: motif for motif in self.motif_list}
        self.motif_name_list = [motif.name for motif in self.motif_list]
        self.n_motifs = len(self.motif_list)
        self.thresholds = {}

        # meta table index is not unique, because one motif may corresponding to multiple gene
        self.meta_table = meta_table

        # add motif cluster and remove duplicates
        if meta_table is not None:
            self.motif_cluster = pd.Series(meta_table[motif_cluster_col].to_dict())
            self.motif_cluster.index.name = "motif"

    def calculate_threshold(self, method="balanced", cpu=1, threshold_value=1000):
        # execute in parallel
        with ProcessPoolExecutor(cpu) as exe:
            futures = {}
            for motif in self.motif_list:
                f = exe.submit(
                    _get_motif_threshold,
                    method=method,
                    motif=motif,
                    threshold_value=threshold_value,
                )
                futures[f] = motif.name

            count = 0
            print("Calculating motif threshold: ", end="")
            for f in as_completed(futures):
                threshold = f.result()
                name = futures[f]
                self.thresholds[name] = threshold
                count += 1
                if count % 250 == 0:
                    print(count, end=" ")
            print()
        return

    def _multi_seq_motif_scores(
            self, seq_list, save_path, region_dim, motif_dim, dtype
    ):
        seq_names = []
        seq_scores = []
        for seq in seq_list:
            this_seq_scores = []
            for motif in self.motif_list:
                threshold = self.thresholds[motif.name]
                n_max_total_score = _single_seq_single_motif_max_score(
                    seq=seq, motif=motif, threshold=threshold
                )
                this_seq_scores.append(n_max_total_score)
            seq_scores.append(this_seq_scores)
            seq_names.append(seq.name)

        seq_scores = np.round(seq_scores)
        max_value = np.iinfo(dtype).max
        seq_scores[seq_scores > max_value] = max_value

        _data = np.array(seq_scores).astype(dtype)
        da = xr.DataArray(
            _data,
            coords=[
                seq_names,
                self.motif_name_list,
                ["n_motifs", "max_score", "total_score"],
            ],
            dims=[region_dim, motif_dim, "motif_value"],
        )
        motif_score_matrix = xr.Dataset({f"{region_dim}_{motif_dim}_da": da})
        motif_score_matrix.to_zarr(save_path)
        return save_path

    def _run_motif_scan_chunks(
            self,
            fasta_path,
            output_dir,
            region_dim,
            motif_dim,
            dtype,
            chunk_size=10000,
            cpu=1,
    ):
        with open(fasta_path) as fasta:
            seqs = []
            for seq in SeqIO.parse(fasta, "fasta"):
                seqs.append(seq)
        print(f"Scan {self.n_motifs} motif in {len(seqs)} sequences.")

        if chunk_size is None:
            chunk_size = min(10000, max(10, len(seqs) // cpu))

        output_dir = pathlib.Path(output_dir).absolute()
        final_dir = f"{output_dir}/{region_dim}_{motif_dim}"
        chunk_dir = f"{output_dir}/.{region_dim}_{motif_dim}_chunk"
        if pathlib.Path(chunk_dir).exists():
            subprocess.run(f"rm -rf {chunk_dir}", shell=True)
        pathlib.Path(chunk_dir).mkdir(parents=True)

        with ProcessPoolExecutor(cpu) as exe:
            futures = {}
            for chunk_i, chunk_start in enumerate(range(0, len(seqs), chunk_size)):
                save_path = f"{chunk_dir}/motif_{chunk_i}.zarr"
                _seqs = seqs[chunk_start: chunk_start + chunk_size]
                future = exe.submit(
                    self._multi_seq_motif_scores,
                    seq_list=_seqs,
                    save_path=str(save_path),
                    region_dim=region_dim,
                    motif_dim=motif_dim,
                    dtype=dtype,
                )
                futures[future] = chunk_i

            chunks_to_write = {}
            for i, future in enumerate(as_completed(futures)):
                print(f"Job {i} returned")
                chunk_i = futures[future]
                output_path = future.result()
                chunks_to_write[chunk_i] = output_path

        write_ordered_chunks(
            chunks_to_write,
            final_path=final_dir,
            append_dim=region_dim,
            engine="zarr",
            coord_dtypes=None,
        )
        subprocess.run(f"rm -rf {chunk_dir}", shell=True)
        return

    def _aggregate_motif_clusters(
            self, output_dir, motif_dim, region_dim, dtype, chunk_size=1000000
    ):
        motif_ds = RegionDS.open(
            f"{output_dir}/{region_dim}_{motif_dim}/", region_dim=region_dim
        )
        _motif_cluster = self.motif_cluster
        _motif_cluster.name = motif_dim
        motif_ds.coords[f"{motif_dim}-cluster"] = _motif_cluster
        motif_cluster_ds = motif_ds.groupby(f"{motif_dim}-cluster").apply(
            lambda _i: _i.max(dim=motif_dim).astype(dtype)
        )

        # save motif cluster matrix by chunk
        motif_cluster_ds = RegionDS(motif_cluster_ds)

        cluster_output_dir = f"{output_dir}/{region_dim}_{motif_dim}-cluster"
        for i, chunk in enumerate(
                motif_cluster_ds.iter_array(
                    dim=region_dim, chunk_size=chunk_size, da="dmr_motif_da", load=True
                )
        ):
            _ds = xr.Dataset({"dmr_motif-cluster_da": chunk}).load()
            if i == 0:
                _ds.to_zarr(cluster_output_dir, mode="w")
            else:
                _ds.to_zarr(cluster_output_dir, append_dim=region_dim)
        return

    def scan_motifs(
            self,
            fasta_path,
            output_dir,
            cpu,
            region_dim,
            motif_dim="motif",
            combine_cluster=True,
            dtype="uint16",
            chunk_size=10000,
    ):
        self._run_motif_scan_chunks(
            fasta_path=fasta_path,
            output_dir=output_dir,
            cpu=cpu,
            region_dim=region_dim,
            motif_dim=motif_dim,
            dtype=dtype,
            chunk_size=chunk_size,
        )
        motif_ds = RegionDS.open(
            f"{output_dir}/{region_dim}_{motif_dim}", region_dim=region_dim
        )

        if combine_cluster:
            self._aggregate_motif_clusters(
                output_dir=output_dir,
                motif_dim=motif_dim,
                region_dim=region_dim,
                dtype=dtype,
                chunk_size=chunk_size,
            )
            motif_cluster_ds = RegionDS.open(
                f"{output_dir}/{region_dim}_{motif_dim}-cluster", region_dim=region_dim
            )
            motif_ds.update(motif_cluster_ds)
        return motif_ds


class Motif(motifs.Motif):
    def __init__(self, alphabet="ACGT", instances=None, counts=None):
        super().__init__(alphabet=alphabet, instances=instances, counts=counts)

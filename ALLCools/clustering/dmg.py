import pathlib

import anndata
import pandas as pd
import xarray as xr
import warnings
import scanpy as sc
import numpy as np
from sklearn.metrics import roc_auc_score
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
from typing import List
from sklearn.metrics import pairwise_distances
from ..mcds import MCDS


def _single_pairwise_dmg(
        cluster_l,
        cluster_r,
        top_n,
        adj_p_cutoff,
        delta_rate_cutoff,
        auroc_cutoff,
        adata_dir,
        dmg_dir,
):
    """Calculate DMG between a pair of adata file"""
    # load data
    adata_l = anndata.read_h5ad(f"{adata_dir}/{cluster_l}.h5ad")
    adata_r = anndata.read_h5ad(f"{adata_dir}/{cluster_r}.h5ad")

    # generate single adata for DMG
    adata = adata_l.concatenate(adata_r, batch_key="groups", index_unique=None)
    adata.obs = pd.concat([adata_l.obs, adata_r.obs])
    try:
        assert adata.obs_names.duplicated().sum() == 0
    except AssertionError as e:
        print(cluster_l, cluster_r)
        raise e

    # reverse_adata, centered by 1 because after normalization all prior is center to 1
    adata.X = (adata.X - 1) * -1 + 1

    # calc DMG
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(
            adata, groupby="groups", n_genes=top_n, method="wilcoxon"
        )
    dmg_result = pd.DataFrame(
        {
            data_key: pd.DataFrame(
                adata.uns["rank_genes_groups"][data_key], columns=[cluster_r, cluster_l]
            ).stack()
            for data_key in ["names", "pvals_adj"]
        }
    )
    dmg_result.reset_index(drop=True, inplace=True)

    # annotate cluster_delta
    dmg_result["left-right"] = f"{cluster_l}-{cluster_r}"
    l_mean = pd.Series(np.mean(adata_l.X, axis=0), index=adata_l.var_names)
    left_mean = dmg_result["names"].map(l_mean)
    r_mean = pd.Series(np.mean(adata_r.X, axis=0), index=adata_l.var_names)
    right_mean = dmg_result["names"].map(r_mean)
    dmg_result["delta"] = left_mean - right_mean

    # filter
    dmg_result = dmg_result[
        (dmg_result["pvals_adj"] < adj_p_cutoff)
        & (dmg_result["delta"].abs() > delta_rate_cutoff)
        ].copy()
    dmg_result["hypo_in"] = dmg_result["delta"].apply(
        lambda i: cluster_l if i < 0 else cluster_r
    )
    dmg_result["hyper_in"] = dmg_result["delta"].apply(
        lambda i: cluster_r if i < 0 else cluster_l
    )
    dmg_result = dmg_result.set_index("names").drop_duplicates()

    # add AUROC and filter again
    auroc = {}
    for gene, row in dmg_result.iterrows():
        yscore = adata.obs_vector(gene)
        ylabel = adata.obs["groups"] == row["hypo_in"]
        score = roc_auc_score(ylabel, yscore)
        score = abs(score - 0.5) + 0.5
        auroc[gene] = score
    dmg_result["AUROC"] = pd.Series(auroc)
    dmg_result = dmg_result[(dmg_result["AUROC"] > auroc_cutoff)].copy()

    # save
    dmg_result.to_hdf(f"{dmg_dir}/{cluster_l}-{cluster_r}.hdf", key="data")
    return


class PairwiseDMG:
    def __init__(
            self,
            max_cell_per_group=1000,
            top_n=10000,
            adj_p_cutoff=0.001,
            delta_rate_cutoff=0.3,
            auroc_cutoff=0.9,
            random_state=0,
            n_jobs=1,
            verbose=True,
    ):
        """
        Calculate pairwise DMGs. After calculation, results saved in self.dmg_table

        Parameters
        ----------
        max_cell_per_group :
            maximum number of cells to use for each group, downsample larger groups to this size
        top_n :
            top N DMGs to report in the final result, if
        adj_p_cutoff :
            adjusted P value cutoff to report significant DMG
        delta_rate_cutoff :
            mC fraction delta cutoff to report significant DMG
        auroc_cutoff :
            AUROC cutoff to report significant DMG
        random_state :
            overall random state to make sure reproducibility
        n_jobs :
            number of cpus
        """
        self.X = None
        self.groups = None
        self._obs_dim = ""
        self._var_dim = ""
        self.dmg_table = None
        self._outlier_label = None
        self.selected_pairs = None

        # parameters
        self.max_cell_per_group = max_cell_per_group
        self.top_n = top_n
        self.adj_p_cutoff = adj_p_cutoff
        self.delta_rate_cutoff = delta_rate_cutoff
        self.auroc_cutoff = auroc_cutoff
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose

        # internal
        self._adata_dir = "_adata_for_dmg"
        self._dmg_dir = "_dmg_results"

    def fit_predict(
            self,
            x,
            groups,
            obs_dim="cell",
            var_dim="gene",
            outlier="Outlier",
            cleanup=True,
            selected_pairs: List[tuple] = None,
    ):
        """
        provide data and perform the pairwise DMG

        Parameters
        ----------
        x :
            2D cell-by-feature xarray.DataArray
        groups :
            cluster labels
        obs_dim :
            name of the cell dim
        var_dim :
            name of the feature dim
        outlier :
            name of the outlier group, if provided, will ignore this label
        cleanup :
            Whether to delete the group adata file
        selected_pairs :
            By default, pairwise DMG will calculate all possible pairs between all the groups, which might be very
            time consuming if the group number is large. With this parameter, you may provide a list of cluster pairs

        Returns
        -------

        """
        if (len(x.shape) != 2) or not isinstance(x, xr.DataArray):
            raise ValueError(
                "Expect an cell-by-feature 2D xr.DataArray as input matrix."
            )
        self._obs_dim = obs_dim
        self._var_dim = var_dim
        self._outlier_label = outlier

        self.X = x
        use_order = x.get_index(obs_dim)
        self.groups = groups.astype("category").reindex(use_order)
        self.selected_pairs = selected_pairs

        # save adata for each group to dict
        if self.verbose:
            print("Generating cluster AnnData files")
        self._save_cluster_adata()

        # run pairwise DMG
        if self.verbose:
            print("Computing pairwise DMG")
        self._pairwise_dmg()

        # cleanup
        if cleanup:
            self._cleanup()

    def _save_cluster_adata(self):
        """Save each group into separate adata, this way reduce the memory during parallel"""
        if self.selected_pairs is not None:
            use_label = set()
            for a, b in self.selected_pairs:
                use_label.add(a)
                use_label.add(b)
        else:
            use_label = None
        adata_dir = pathlib.Path(self._adata_dir)
        adata_dir.mkdir(exist_ok=True)

        total_cells = self.X.get_index(self._obs_dim)

        for cluster, sub_series in self.groups.groupby(self.groups):
            if (use_label is not None) and (cluster not in use_label):
                continue
            if cluster == self._outlier_label:
                # skip outlier
                continue
            output_path = adata_dir / f"{cluster}.h5ad"
            if output_path.exists():
                continue
            sub_series = sub_series.cat.remove_unused_categories()
            if sub_series.size > self.max_cell_per_group:
                sub_series = sub_series.sample(
                    self.max_cell_per_group, random_state=self.random_state
                )
            cluster_adata = anndata.AnnData(
                X=self.X.sel(
                    {self._obs_dim: total_cells.intersection(sub_series.index)}
                ).values,
                obs=pd.DataFrame({"groups": sub_series.astype("category")}),
                var=pd.DataFrame([], index=self.X.get_index(self._var_dim)),
            )
            cluster_adata.write_h5ad(output_path)
        return

    def _pairwise_dmg(self):
        """pairwise DMG runner, result save to self.dmg_table"""
        dmg_dir = pathlib.Path(self._dmg_dir)
        dmg_dir.mkdir(exist_ok=True)

        if self.selected_pairs is None:
            pairs = [
                i
                for i in combinations(sorted(self.groups.unique()), 2)
                if self._outlier_label not in i
            ]
        else:
            if self.verbose:
                print("Using user provided gene pairs")
            pairs = self.selected_pairs
        if self.verbose:
            print(len(pairs), "pairwise DMGs")
        n_pairs = len(pairs)

        with ProcessPoolExecutor(min(n_pairs, self.n_jobs)) as exe:
            step = max(1, n_pairs // 10)
            futures = {}
            for (cluster_l, cluster_r) in pairs:
                f = exe.submit(
                    _single_pairwise_dmg,
                    cluster_l=cluster_l,
                    cluster_r=cluster_r,
                    top_n=self.top_n,
                    adj_p_cutoff=self.adj_p_cutoff,
                    delta_rate_cutoff=self.delta_rate_cutoff,
                    auroc_cutoff=self.auroc_cutoff,
                    adata_dir=self._adata_dir,
                    dmg_dir=self._dmg_dir,
                )
                futures[f] = (cluster_l, cluster_r)
            for i, f in enumerate(as_completed(futures)):
                f.result()
                if i % step == 0 and self.verbose:
                    print(f"{i + 1}/{n_pairs} finished")

        # summarize
        self.dmg_table = pd.concat((pd.read_hdf(p) for p in dmg_dir.glob("*.hdf")))
        return

    def _cleanup(self):
        """Delete group adata files"""
        subprocess.run(
            ["rm", "-rf", str(self._adata_dir), str(self._dmg_dir)], check=True
        )

    def aggregate_pairwise_dmg(self, adata, groupby, obsm="X_pca"):
        """
        Aggregate pairwise DMG results for each cluster, rank DMG for the cluster by the sum of
        AUROC * cluster_pair_similarity
        This way, the DMGs having large AUROC between similar clusters get more weights

        Parameters
        ----------
        adata :
        groupby :
        obsm :

        Returns
        -------

        """
        # using the cluster centroids in PC space to calculate dendrogram
        pc_matrix = adata.obsm[obsm]
        pc_center = (
            pd.DataFrame(pc_matrix, index=adata.obs_names)
                .groupby(adata.obs[groupby])
                .median()
        )
        # calculate cluster pairwise similarity
        cluster_dist = pairwise_distances(pc_center)
        cluster_dist = pd.DataFrame(
            cluster_dist, index=pc_center.index, columns=pc_center.index
        )
        cluster_dist_norm = cluster_dist / cluster_dist.values.max()
        cluster_sim = 1 - cluster_dist_norm
        cluster_pair_sim_dict = {
            f"{a}-{b}": value for (a, b), value in cluster_sim.unstack().iteritems()
        }
        self.dmg_table["similarity"] = self.dmg_table["left-right"].map(
            cluster_pair_sim_dict
        )

        # aggregate pairwise DMG to get the cluster level DMG, use the similarity to normalize AUROC
        cluster_dmgs = {}
        for cluster, sub_df in self.dmg_table.groupby("hypo_in"):
            values = sub_df["AUROC"] * sub_df["similarity"]
            cluster_dmg_order = (
                values.groupby(values.index).sum().sort_values(ascending=False)
            )
            cluster_dmgs[cluster] = cluster_dmg_order
        return cluster_dmgs


def _single_ovr_dmg(
        cell_label,
        mcds,
        obs_dim,
        var_dim,
        mc_type,
        top_n,
        adj_p_cutoff,
        fc_cutoff,
        auroc_cutoff,
):
    """single one vs rest DMG runner"""
    # get adata
    cell_judge = mcds.get_index(obs_dim).isin(cell_label.index)
    adata = mcds.sel({obs_dim: cell_judge}).get_adata(
        mc_type=mc_type, var_dim=var_dim, select_hvf=False
    )
    adata.var = pd.DataFrame([], index=adata.var_names)
    adata.obs["groups"] = cell_label.astype("category")

    # calc DMG
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(
            adata, groupby="groups", n_genes=top_n, method="wilcoxon"
        )
    dmg_result = pd.DataFrame(
        {
            data_key: pd.DataFrame(adata.uns["rank_genes_groups"][data_key]).stack()
            for data_key in ["names", "pvals_adj"]
        }
    )
    dmg_result = dmg_result[
        dmg_result.index.get_level_values(1).astype(bool)
    ].reset_index(drop=True)
    # add fold change
    in_cells_mean = adata.X[
        adata.obs["groups"].astype(bool),
    ].mean(axis=0)
    out_cells_mean = adata.X[
        ~adata.obs["groups"].astype(bool),
    ].mean(axis=0)
    fc = pd.Series(in_cells_mean / out_cells_mean, index=adata.var_names)
    dmg_result["fc"] = dmg_result["names"].map(fc)
    # filter
    dmg_result = dmg_result[
        (dmg_result["pvals_adj"] < adj_p_cutoff) & (dmg_result["fc"] < fc_cutoff)
        ].copy()
    dmg_result = dmg_result.set_index("names").drop_duplicates()

    # add AUROC and filter again
    auroc = {}
    for gene, row in dmg_result.iterrows():
        yscore = adata.obs_vector(gene)
        ylabel = adata.obs["groups"] == True
        score = roc_auc_score(ylabel, yscore)
        score = abs(score - 0.5) + 0.5
        auroc[gene] = score
    dmg_result["AUROC"] = pd.Series(auroc)
    dmg_result = dmg_result[(dmg_result["AUROC"] > auroc_cutoff)].copy()
    return dmg_result


def _one_vs_rest_dmr_runner(
        cell_meta,
        group,
        cluster,
        max_cluster_cells,
        max_other_fold,
        mcds_paths,
        obs_dim,
        var_dim,
        mc_type,
        top_n,
        adj_p_cutoff,
        fc_cutoff,
        auroc_cutoff,
        verbose = True,
):
    """one vs rest DMG runner"""
    if verbose:
        print(f"Calculating cluster {cluster} DMGs.")

    mcds = MCDS.open(mcds_paths)
    # determine cells to use
    cluster_judge = cell_meta[group] == cluster
    in_cells = cluster_judge[cluster_judge]
    out_cells = cluster_judge[~cluster_judge]
    if in_cells.size > max_cluster_cells:
        in_cells = in_cells.sample(max_cluster_cells, random_state=0)

    max_other_cells = in_cells.size * max_other_fold
    if out_cells.size > max_other_cells:
        out_cells = out_cells.sample(max_other_cells, random_state=0)

    cell_label = pd.concat([in_cells, out_cells])
    dmg_df = _single_ovr_dmg(
        cell_label=cell_label,
        mcds=mcds,
        obs_dim=obs_dim,
        var_dim=var_dim,
        mc_type=mc_type,
        top_n=top_n,
        adj_p_cutoff=adj_p_cutoff,
        fc_cutoff=fc_cutoff,
        auroc_cutoff=auroc_cutoff,
    )
    return dmg_df


def one_vs_rest_dmg(
        cell_meta,
        group,
        mcds=None,
        mcds_paths=None,
        obs_dim="cell",
        var_dim="gene",
        mc_type="CHN",
        top_n=1000,
        adj_p_cutoff=0.01,
        fc_cutoff=0.8,
        auroc_cutoff=0.8,
        max_cluster_cells=2000,
        max_other_fold=5,
        cpu=1,
        verbose=True,
):
    """
    Calculating cluster marker genes using one-vs-rest strategy.

    Parameters
    ----------
    cell_meta
        cell metadata containing cluster labels
    group
        the name of the cluster label column
    mcds
        cell-by-gene MCDS object for calculating DMG. Provide either mcds_paths or mcds.
    mcds_paths
        cell-by-gene MCDS paths for calculating DMG. Provide either mcds_paths or mcds.
    obs_dim
        dimension name of the cells
    var_dim
        dimension name of the features
    mc_type
        value to select methylation type in the mc_type dimension
    top_n
        report top N DMGs
    adj_p_cutoff
        adjusted P value cutoff to report significant DMG
    fc_cutoff
        mC fraction fold change cutoff to report significant DMG
    auroc_cutoff
        AUROC cutoff to report significant DMG
    max_cluster_cells
        The maximum number of cells from a group, downsample large group to this number
    max_other_fold
        The fold of other cell numbers comparing
    cpu :
            number of cpus
    Returns
    -------
    dmg_table
        pandas Dataframe of the one-vs-rest DMGs
    """
    tmp = None
    if mcds_paths is not None:
        tmp_created = False
    elif mcds is not None:
        tmp = 'tmp_one_vs_rest.mcds'
        mcds.to_zarr(tmp)
        mcds_paths = tmp
        tmp_created = True
    else:
        raise ValueError(f'Need to provide either "mcds_path" or "mcds".')

    clusters = cell_meta[group].unique()
    dmg_table = []
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for cluster in sorted(clusters):
            f = exe.submit(
                _one_vs_rest_dmr_runner,
                cell_meta=cell_meta,
                group=group,
                cluster=cluster,
                max_cluster_cells=max_cluster_cells,
                max_other_fold=max_other_fold,
                mcds_paths=mcds_paths,
                obs_dim=obs_dim,
                var_dim=var_dim,
                mc_type=mc_type,
                top_n=top_n,
                adj_p_cutoff=adj_p_cutoff,
                fc_cutoff=fc_cutoff,
                auroc_cutoff=auroc_cutoff,
            )
            futures[f] = cluster

        for f in as_completed(futures):
            cluster = futures[f]
            if verbose:
                print(f"{cluster} Finished.")
            dmg_df = f.result()
            dmg_df["cluster"] = cluster
            dmg_table.append(dmg_df)
    dmg_table = pd.concat(dmg_table)

    if tmp_created:
        import shutil
        shutil.rmtree(tmp)

    return dmg_table

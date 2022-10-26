import numpy as np
import pandas as pd
from sklearn.preprocessing import scale


def _mvalue(da, scale_group=True, alpha=1):
    mc = da.sel(count_type="mc")
    cov = da.sel(count_type="cov")
    m = np.log2((mc + alpha) / (cov - mc + alpha))
    if scale_group:
        group_axis = da.dims.index("group")
        m.values = scale(m, axis=group_axis)
    return m


def _frac(da, alpha=0.01):
    f = da.sel(count_type="mc") / (da.sel(count_type="cov") + alpha)
    return f


class DMSAggregate:
    def __init__(
        self,
        base_ds,
        cell_groups,
        dms_ds,
        chrom,
        start,
        end,
        mcg_pattern="CGN",
        sample_dim="sample_id",
        count_da="data",
        quantile=0.1,
        min_delta=0.2,
        frac_delta_cutoff=0.2,
        frac_delta_quantile=0.05,
        max_dist=250,
        corr_cutoff=0.3,
        drop_additional_samples=False,
    ):
        self.base_ds = base_ds
        self.cell_groups = cell_groups
        self.dms_ds = dms_ds
        self.chrom = chrom
        self.start = start
        self.end = end
        self.mcg_pattern = mcg_pattern
        self.sample_dim = sample_dim
        self.count_da = count_da

        self.quantile = quantile
        self.min_delta = min_delta
        self.frac_delta_cutoff = frac_delta_cutoff
        self.frac_delta_quantile = frac_delta_quantile
        self.max_dist = max_dist
        self.corr_cutoff = corr_cutoff
        self.drop_additional_samples = drop_additional_samples
        self._all_cpgs = None
        return

    def _get_all_cpg(self):
        if self.drop_additional_samples:
            use_samples = pd.Series(self.cell_groups).index
            sample_index = self.base_ds.get_index(self.sample_dim)
            self.base_ds = self.base_ds.sel({self.sample_dim: sample_index.isin(use_samples)})

        all_cpgs = self.base_ds.fetch(self.chrom, self.start, self.end).select_mc_type(self.mcg_pattern)

        all_cpgs.coords["group"] = all_cpgs[self.sample_dim].to_pandas().map(self.cell_groups)
        all_cpgs["group_data"] = (
            all_cpgs[self.count_da].groupby("group").sum(dim=self.sample_dim).load(scheduler="sync")
        )

        self._all_cpgs = all_cpgs
        return

    def _add_group_value_to_dms_ds(self):
        dms_ds = self.dms_ds

        # TODO check dms_ds pos order, make sure its always ordered
        sort_pos = dms_ds.get_index("pos").sort_values()
        dms_ds = dms_ds.sel(pos=sort_pos)

        dms_count = self._all_cpgs.fetch_positions(positions=dms_ds.get_index("pos")).load(scheduler="sync")
        dms_ds["group_data"] = dms_count["group_data"]

        # get both forward and reverse strand
        dms_pos = pd.Index(dms_ds.get_index("pos").tolist() + (dms_ds.get_index("pos") + 1).tolist()).sort_values()
        self._all_cpgs.coords["is_dms"] = self._all_cpgs.get_index("pos").isin(dms_pos)

        dms_ds["mvalue"] = _mvalue(dms_ds["group_data"])
        dms_ds["frac"] = _frac(dms_ds["group_data"])

        self.dms_ds = dms_ds
        return

    def _filter_by_frac_delta(self):
        _q = min(self.frac_delta_quantile, 1 - self.frac_delta_quantile)
        frac_delta = self.dms_ds["frac"].quantile((1 - _q), dim="group") - self.dms_ds["frac"].quantile(_q, dim="group")
        frac_delta_judge = frac_delta > self.frac_delta_cutoff
        self.dms_ds = self.dms_ds.sel(pos=frac_delta_judge)

    def _combine_dms_to_dmr(self):
        dms_ds = self.dms_ds

        # step 1: combine raw DMR windows based on distance
        dms_dist = dms_ds["pos"][1:].values - dms_ds["pos"][:-1].values
        dist_judge = dms_dist < self.max_dist
        cur_dmr = 0
        dmr_ids = [0]
        for i in dist_judge:
            if not i:
                cur_dmr += 1
            dmr_ids.append(cur_dmr)
        dmr_ids = pd.Series(dmr_ids, index=dms_ds.get_index("pos"))

        # step 2: determine correlation between adjacent DMS
        data = dms_ds["mvalue"].transpose("pos", "group").to_pandas()
        a = data.iloc[:-1, :].reset_index(drop=True)
        b = data.iloc[1:, :].reset_index(drop=True)
        # index of corr means the correlation of that DMS with the previous DMS
        # regardless of the genome distance
        # fill na value with 1, tend to merge DMR due to nan
        corrs = a.fillna(1).corrwith(b.fillna(1), axis=1, method="pearson")
        corrs.index = data.iloc[1:, :].index.copy()

        # step 3: recalculate DMR windows based on both distance and correlation
        raw_dmr_table = pd.DataFrame({"dmr": dmr_ids, "corr": corrs})
        dmr_dict = {}
        cur_dmr_id = 0
        for _, sub_df in raw_dmr_table.groupby("dmr"):
            dmr_dict[sub_df.index[0]] = cur_dmr_id
            if sub_df.shape[0] > 1:
                for dms_id, corr in sub_df["corr"][1:].items():
                    if corr > self.corr_cutoff:
                        dmr_dict[dms_id] = cur_dmr_id
                    else:
                        cur_dmr_id += 1
                        dmr_dict[dms_id] = cur_dmr_id
            cur_dmr_id += 1
        dmr_ids = pd.Series(dmr_dict)
        dmr_ids.index.name = "pos"
        dms_ds.coords["dmr"] = dmr_ids

        self.dms_ds = dms_ds
        return

    def _get_dmr_ds(self):
        # extend DMR ids to both strand
        dmr_ids = self.dms_ds.coords["dmr"].to_pandas()
        dmr_ids_reverse = dmr_ids.copy()
        dmr_ids_reverse.index += 1
        dmr_ids = pd.concat([dmr_ids, dmr_ids_reverse]).sort_index()

        dms_sample_count = self._all_cpgs[["data"]].sel(pos=dmr_ids.index)
        dms_sample_count.coords["dmr"] = dmr_ids
        dmr_sample_count = dms_sample_count.groupby("dmr").sum(dim="pos")

        # get DMR region
        dmr_region = (
            dmr_ids.groupby(dmr_ids)
            .apply(lambda i: pd.Series({"start": i.index.min(), "end": i.index.max() + 1}))
            .unstack()
        )
        dmr_region.index.name = "dmr"
        dmr_sample_count.coords["dmr_start"] = dmr_region["start"]
        dmr_sample_count.coords["dmr_end"] = dmr_region["end"]
        n_dmr = dmr_sample_count.get_index("dmr").size
        dmr_sample_count.coords["dmr_chrom"] = pd.Series(
            [self._all_cpgs.chrom] * n_dmr, index=dmr_sample_count.get_index("dmr")
        )
        dmr_sample_count["data"] = dmr_sample_count["data"].astype("uint32").load(scheduler="sync")
        dmr_sample_count.coords["dmr_length"] = dmr_region["end"] - dmr_region["start"]

        dmr_sample_count["dmr"] = (
            dmr_sample_count["dmr_chrom"].to_pandas()
            + "-"
            + dmr_sample_count["dmr_start"].astype(str).to_pandas()
            + "-"
            + dmr_sample_count["dmr_length"].astype(str).to_pandas()
        )
        return dmr_sample_count

    def run(self):
        self._get_all_cpg()
        self._add_group_value_to_dms_ds()
        self._filter_by_frac_delta()
        self._combine_dms_to_dmr()
        dmr_ds = self._get_dmr_ds()
        return dmr_ds

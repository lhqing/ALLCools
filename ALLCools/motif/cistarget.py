"""
Reimplement the core cisTarget function for motif enrichment analysis

Code is adapted from pycistarget package, the ctxcore package contains core functions for cisTarget

pycistarget LICENSE
https://github.com/aertslab/pycistarget/blob/master/LICENSE.txt
"""

import numpy as np
import pandas as pd
from ctxcore.recovery import aucs as calc_aucs
from ctxcore.recovery import recovery


def cistarget_motif_enrichment(
    rank_df, total_regions, auc_threshold=0.005, nes_threshold=3, rank_threshold=0.05, full_motif_stats=False
):
    """
    Perform the cistarget motif enrichment analysis by considering normalized motif score rank AUC.

    Parameters
    ----------
    rank_df :
        The motif score rank dataframe.
    total_regions :
        The total number of regions.
    auc_threshold :
        The total proportion of regions to be considered when calculating AUC.
    nes_threshold :
        The normalized enrichment score threshold to determine if a motif is enriched.
    rank_threshold :
        The rank threshold to determine if a motif hit in a region.
    full_motif_stats :
        Whether to return the full motif statistics.

    Returns
    -------
    motif_enrichment :
        The motif enrichment dataframe.
    motif_hits :
        The motif hits in regions.
    """
    # Get features, rankings and weights
    motifs, _ = rank_df.index.values, rank_df.values

    weights = np.asarray(np.ones(rank_df.shape[1]))

    # Calculate recovery curves, AUC and NES values.
    aucs = calc_aucs(
        rnk=rank_df,
        total_genes=total_regions,
        weights=weights,
        auc_threshold=auc_threshold,
    )
    ness = (aucs - aucs.mean()) / aucs.std()

    # Keep only features that are enriched, i.e. NES sufficiently high.
    enriched_features_idx = ness >= nes_threshold
    weights = np.asarray(np.ones(rank_df.shape[1]))

    # terminate if no features are enriched
    # noinspection PyTypeChecker
    if sum(enriched_features_idx) == 0:
        print("No enriched motifs found")
        motif_enrichment = pd.DataFrame([], columns=["NES", "AUC", "rank_at_max", "n_motif_hit"]).astype(
            {"NES": float, "AUC": float, "rank_at_max": int, "n_motif_hit": int}
        )
        motif_hits = pd.DataFrame([], columns=rank_df.columns, dtype=bool)
        return motif_enrichment, motif_hits

    # Make dataframe
    enriched_features = pd.DataFrame(
        index=pd.Index(motifs[enriched_features_idx], name="MotifID"),
        data={
            "NES": ness[enriched_features_idx],
            "AUC": aucs[enriched_features_idx],
        },
    )
    # full motif stats
    if full_motif_stats:
        full_stats = pd.DataFrame(
            index=pd.Index(motifs, name="MotifID"),
            data={
                "NES": ness,
                "AUC": aucs,
            },
        )
    else:
        full_stats = None

    # Recovery analysis (Leading edge analysis)
    rccs, _ = recovery(
        rnk=rank_df,
        total_genes=total_regions,
        weights=weights,
        rank_threshold=int(rank_threshold * total_regions),
        auc_threshold=auc_threshold,
        no_auc=True,
    )
    avg_rcc = rccs.mean(axis=0)
    avg_2std_rcc = avg_rcc + 2.0 * rccs.std(axis=0)
    # Select max rank for each motif, smaller rank number is considered motif hit
    rccs = rccs[enriched_features_idx, :]
    rank_max_record = {}
    for i, feature in enumerate(enriched_features.index):
        rcc = rccs[i]
        rank_at_max = np.argmax(rcc - avg_2std_rcc)
        rank_max_record[feature] = rank_at_max

    # Format enriched features
    enriched_features["rank_at_max"] = pd.Series(rank_max_record)
    enriched_features = enriched_features.sort_values("NES", ascending=False)
    motif_enrichment = enriched_features[["NES", "AUC", "rank_at_max"]]

    motif_hit = rank_df.loc[motif_enrichment.index] < motif_enrichment["rank_at_max"].values[:, None]
    # noinspection PyUnresolvedReferences
    motif_enrichment["n_motif_hit"] = motif_hit.sum(axis=1)
    return motif_enrichment, motif_hit, full_stats

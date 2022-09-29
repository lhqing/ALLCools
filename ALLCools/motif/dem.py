import numpy as np
import pandas as pd
from scipy.stats import ranksums
from sklearn.metrics import roc_curve
from statsmodels.stats.multitest import multipletests


def _get_optimal_threshold(scores, labels):
    pos_scores = scores > 0
    _labels = labels[pos_scores]
    _scores = scores[pos_scores]
    fpr, tpr, thresholds = roc_curve(_labels, _scores)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    return optimal_threshold


def dem_motif_enrichment(hypo_score_df, hyper_score_df, log2_fc_threshold=0.5, alpha=0.05, motif_score_threshold=3):
    """
    Perform the Wilcoxon rank sum test on hypo and hyper region motif scores.

    Parameters
    ----------
    hypo_score_df :
        The motif score dataframe of hypo regions.
    hyper_score_df :
        The motif score dataframe of hyper regions.
    log2_fc_threshold :
        The log2 fold change threshold to determine if a motif is enriched.
    alpha :
        The adjusted p-value threshold to determine if a motif is enriched.
    motif_score_threshold :
        The motif score threshold to determine if a motif hit in a region.

    Returns
    -------
    motif_enrichment :
        The motif enrichment dataframe.
    hypo_motif_hits :
        The motif hits in hypo regions.
    hyper_motif_hits :
        The motif hits in hyper regions.
    """
    # determine hypo hyper regions are not overlapped
    assert hypo_score_df.columns.intersection(hyper_score_df.columns).size == 0, "hypo and hyper regions are overlapped"

    # determine motif rows are the same
    # noinspection PyTypeChecker
    assert sum(hypo_score_df.index != hyper_score_df.index) == 0, "motif rows are not the same"

    fg_mat = hypo_score_df.values
    bg_mat = hyper_score_df.values
    motifs = hypo_score_df.index
    wilcox_test = [ranksums(fg_mat[x], y=bg_mat[x]) for x in range(fg_mat.shape[0])]

    # Log2FC
    mean_fg = fg_mat.mean(axis=1)
    mean_bg = bg_mat.mean(axis=1)
    log_fc = np.log2((mean_fg + 10**-12) / (mean_bg + 10**-12)).tolist()

    # P-value correction
    # noinspection PyUnresolvedReferences
    p_value = [w.pvalue for w in wilcox_test]

    judge, q_value, *_ = multipletests(p_value, alpha=alpha)

    # Motif df
    motif_df = pd.DataFrame(
        {"log2_fc": log_fc, "q_value": q_value, "mean_fg": mean_fg, "mean_bg": mean_bg}, index=motifs
    )

    fc_judge = motif_df["log2_fc"].abs() > log2_fc_threshold
    motif_df = motif_df[fc_judge & judge].copy()

    # Motif hits versus background
    keep_motifs = motif_df.index.tolist()
    keep_motifs_bool = motifs.isin(keep_motifs)

    scores_mat = np.concatenate([fg_mat[keep_motifs_bool], bg_mat[keep_motifs_bool]], axis=1)
    labels = np.repeat([1, 0], (fg_mat.shape[1], bg_mat.shape[1]))
    motif_hit_thresholds = []
    for i in range(scores_mat.shape[0]):
        opt_score = _get_optimal_threshold(scores=scores_mat[i], labels=labels)
        hit_threshold = max(motif_score_threshold, opt_score)
        motif_hit_thresholds.append(hit_threshold)
    motif_hit_thresholds = np.array(motif_hit_thresholds)

    hypo_motif_hits = pd.DataFrame(
        fg_mat[keep_motifs_bool] > motif_hit_thresholds[:, None], index=keep_motifs, columns=hypo_score_df.columns
    )
    hyper_motif_hits = pd.DataFrame(
        bg_mat[keep_motifs_bool] > motif_hit_thresholds[:, None], index=keep_motifs, columns=hyper_score_df.columns
    )
    return motif_df, hypo_motif_hits, hyper_motif_hits

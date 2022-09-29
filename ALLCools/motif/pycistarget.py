# """
# pycistarget LICENSE
# https://github.com/aertslab/pycistarget/blob/master/LICENSE.txt
# """
# from functools import partial
# from itertools import repeat
#
# import numpy as np
# import pandas as pd
# from ctxcore.recovery import aucs as calc_aucs
# from ctxcore.recovery import leading_edge4row, recovery
#
#
# def cistarget_motif_enrichment(db_rankings_regions, total_regions, auc_threshold, nes_threshold, rank_threshold, name):
#     # Get features, rankings and weights
#     features, rankings = db_rankings_regions.index.values, db_rankings_regions.values
#
#     weights = np.asarray(np.ones(db_rankings_regions.shape[1]))
#     # Calculate recovery curves, AUC and NES values.
#     aucs = calc_aucs(db_rankings_regions, total_regions, weights, auc_threshold)
#     ness = (aucs - aucs.mean()) / aucs.std()
#     # Keep only features that are enriched, i.e. NES sufficiently high.
#     enriched_features_idx = ness >= nes_threshold
#     # terminate if no features are enriched
#     # noinspection PyTypeChecker
#     if sum(enriched_features_idx) == 0:
#         print(f"No enriched motifs found for {name}")
#         motif_enrichment = pd.DataFrame(
#             data={
#                 "Logo": [],
#                 "Region_set": [],
#                 "Direct_annot": [],
#                 "Motif_similarity_annot": [],
#                 "Orthology_annot": [],
#                 "Motif_similarity_and_Orthology_annot": [],
#                 "NES": [],
#                 "AUC": [],
#                 "Rank_at_max": [],
#             }
#         )
#         motif_hits = {"Database": {}, "Region_set": {}}
#         return motif_enrichment, motif_hits
#
#     COLUMN_NAME_NES = "NES"
#     COLUMN_NAME_AUC = "AUC"
#     COLUMN_NAME_GRP = "GROUP"
#     COLUMN_NAME_MOTIF_ID = "MotifID"
#     COLUMN_NAME_TARGET_GENES = "TargetRegions"
#     COLUMN_NAME_RANK_AT_MAX = "RankAtMax"
#
#     # Make dataframe
#     enriched_features = pd.DataFrame(
#         index=pd.Index(features[enriched_features_idx], name=COLUMN_NAME_MOTIF_ID),
#         data={
#             COLUMN_NAME_NES: ness[enriched_features_idx],
#             COLUMN_NAME_AUC: aucs[enriched_features_idx],
#             COLUMN_NAME_GRP: repeat(region_set_signature.transcription_factor, sum(enriched_features_idx)),
#         },
#     )
#
#     # Recovery analysis
#     rccs, _ = recovery(
#         db_rankings_regions, total_regions, weights, int(rank_threshold * total_regions), auc_threshold, no_auc=True
#     )
#     avgrcc = rccs.mean(axis=0)
#     avg2stdrcc = avgrcc + 2.0 * rccs.std(axis=0)
#
#     # Select features
#     rccs = rccs[enriched_features_idx, :]
#     rankings = rankings[enriched_features_idx, :]
#
#     # Format df
#     enriched_features.columns = pd.MultiIndex.from_tuples(list(zip(repeat("Enrichment"), enriched_features.columns)))
#     df_rnks = pd.DataFrame(
#         index=enriched_features.index,
#         columns=pd.MultiIndex.from_tuples(list(zip(repeat("Ranking"), regions))),
#         data=rankings,
#     )
#     df_rccs = pd.DataFrame(
#         index=enriched_features.index,
#         columns=pd.MultiIndex.from_tuples(
#             list(zip(repeat("Recovery"), np.arange(int(rank_threshold * total_regions))))
#         ),
#         data=rccs,
#     )
#     enriched_features = pd.concat([enriched_features, df_rccs, df_rnks], axis=1)
#
#     # Calculate the leading edges for each row. Always return importance from gene inference phase.
#     weights = np.asarray([region_set_signature[region] for region in regions])
#     enriched_features[
#         [("Enrichment", COLUMN_NAME_TARGET_GENES), ("Enrichment", COLUMN_NAME_RANK_AT_MAX)]
#     ] = enriched_features.apply(
#         partial(leading_edge4row, avg2stdrcc=avg2stdrcc, genes=regions, weights=weights), axis=1
#     )
#     enriched_features = enriched_features["Enrichment"].rename_axis(None)
#
#     # Format enriched features
#     enriched_features.columns = ["NES", "AUC", "Region_set", "Motif_hits", "Rank_at_max"]
#     enriched_features = enriched_features.sort_values("NES", ascending=False)
#     motif_enrichment = enriched_features[["Region_set", "NES", "AUC", "Rank_at_max"]]
#
#     # Motif hits
#     db_motif_hits = {
#         key: [
#             enriched_features.loc[key, "Motif_hits"][i][0] for i in range(len(enriched_features.loc[key, "Motif_hits"]))
#         ]
#         for key in enriched_features.index
#     }
#     motif_hits = {"Database": db_motif_hits}
#     motif_enrichment["Motif_hits"] = [len(db_motif_hits[i]) for i in db_motif_hits.keys()]
#     return

from collections import Counter

import numpy as np
import pandas as pd


def mc_gene_panel_greedy_search(cluster_cef_path,
                                pwdmg_path,
                                cef_value_ascending,
                                n_cef=500,
                                final_genes=None,
                                min_dmg_support=2,
                                selection_target=100,
                                add_gene_step=10):
    if final_genes is None:
        final_genes = {}
    else:
        # genes are preselected
        if not isinstance(final_genes, dict):
            final_genes = {g: 'pre-selected' for g in final_genes}

    # only select top n_cef genes for each cluster, when the cluster is involved in PWDMG
    cef = pd.read_csv(cluster_cef_path, index_col=0)
    selected_cef = {col: cef[col].dropna().sort_values(ascending=cef_value_ascending)[:n_cef].index
                    for col in cef.columns}

    # for each cluster pairs, save all significant DMGs for selection
    pwdmg = pd.read_csv(pwdmg_path, index_col=0)
    selected_genes = {}
    for (a, b), sub_df in pwdmg.groupby(['hypo_in', 'hyper_in']):
        possible_gene = sub_df.loc[sub_df.index.isin(selected_cef[a])].index
        selected_genes[(a, b)] = possible_gene

    # greedy search genes satisfy most PWDMG pairs, until any of these conditions achieved
    # 1. all pairs have > min_dmg_support genes included
    # 2. the selection_target gene number reached
    # 3. there is no more significant PWDMG & CEF remained
    final_genes, unsatisfied_pairs = dict_greedy_search(
        final_dict=final_genes,
        source_dict=selected_genes,
        min_value_support=min_dmg_support,
        selection_target=selection_target,
        add_gene_step=add_gene_step)

    # print final info after selection
    print(f'{len(final_genes)} genes in final list.')
    if len(unsatisfied_pairs) == 0:
        print(f'All cluster pairs satisfied by having at least {min_dmg_support} DMGs')
    else:
        print(f'{len(unsatisfied_pairs)} cluster pairs unsatisfied.')
        print(unsatisfied_pairs)
        for c, i in pd.Series(np.concatenate(unsatisfied_pairs)).value_counts().items():
            print(f'Cluster {c} in {i} pairs')
    return final_genes


def dict_greedy_search(source_dict,
                       final_dict: dict,
                       min_value_support,
                       selection_target,
                       add_gene_step):
    dmg_rank = 0
    while True:
        unsatisfied_keys = []
        cur_vote = Counter()
        for key, genes in source_dict.items():
            if genes.isin(final_dict).sum() >= min_value_support:
                continue
            else:
                cur_vote += Counter(genes)
                unsatisfied_keys.append(key)
        if (len(cur_vote) > 0) and (len(final_dict) < selection_target):
            _this_step = min(add_gene_step, selection_target - len(final_dict))
            ordered_genes = pd.Series(cur_vote).sort_values(ascending=False).index
            genes_to_add = 0
            reasons = {}
            for g in ordered_genes:
                if g in final_dict:
                    continue
                genes_to_add += 1
                if genes_to_add > _this_step:
                    break
                reasons[g] = f'Greedy-{dmg_rank}'
                dmg_rank += 1
            if genes_to_add < _this_step:
                # no enough genes to add
                print('no enough genes to add')
                break
            final_dict.update(reasons)
        else:
            break
    return final_dict, unsatisfied_keys

from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from networkx.algorithms.centrality import edge_betweenness_centrality
from numpy import log
from scipy.special import betaln

from .dendrogram import extract_all_nodes
from ALLCools.plot.dendro import *


def linkage_to_graph(linkage):
    """Turn the linkage matrix into a graph, an epimutation will just be remove one edge from the graph"""
    _linkage = linkage.astype(int)
    n_leaf = _linkage.shape[0] + 1
    edges = []
    for i in range(_linkage.shape[0]):
        cur_node = i + n_leaf
        left, right, *_ = _linkage.iloc[i]
        edges.append([left, cur_node])
        edges.append([right, cur_node])
    g = nx.Graph()
    g.add_edges_from(edges)
    return g


def cut_by_highest_betweenness_centrality(g):
    # order graph node by betweenness_centrality
    highest_centrality_edge = pd.Series(edge_betweenness_centrality(g)).sort_values(ascending=False).index[0]
    _g = g.copy()
    _g.remove_edge(*highest_centrality_edge)
    left_tree, right_tree = nx.connected_component_subgraphs(_g)
    return left_tree, right_tree, highest_centrality_edge


def log_proba_beta_binomial(x, n, a, b):
    """log likelihood for the beta-binomial dist, ignore part not related to a and b."""
    like = betaln((a + x), (b + n - x)) - betaln(a, b)
    # when a or b has 0, like will have nan
    return like.fillna(0)


def parse_one_pattern(tree_g, edges_to_remove, mc_df, cov_df):
    """
    for a particular epimutation combination (edges_to_remove),
    calculate the a and b for beta-binomial dist in each leaf node group.

    after removing the edges (epimutations),
    the leaf node group are leaf nodes in each of the disconnected sub graph.
    """
    group_mc_df = mc_df.copy()
    group_un_mc_df = cov_df - group_mc_df

    sub_g = tree_g.copy()
    if len(edges_to_remove) > 0:  # this is the case of adding empty edge by left-right combine
        sub_g.remove_edges_from(edges_to_remove)
    # get disconnected sub-graphs
    sub_tree = nx.connected_component_subgraphs(sub_g)

    # for each sub-graph, add up the mc and un-mc of all leaf nodes for group a, b in beta-binomial dist
    for _tree in sub_tree:
        judge = group_mc_df.columns.isin(_tree.nodes)
        if judge.sum() == 0:
            # if sub-graph do not have leaf nodes, skip this sub-graph
            continue
        group_mc_df.loc[:, judge] = group_mc_df.loc[:, judge].sum(
            axis=1).values[:, None]
        group_un_mc_df.loc[:, judge] = group_un_mc_df.loc[:, judge].sum(
            axis=1).values[:, None]

    # group_mc_df is a, group_un_mc_df is b for beta-binomial dist
    # each group of leaf nodes share same a, b
    return group_mc_df, group_un_mc_df


def mutation_likelihood(n_mutation, p_mutation, n_edges):
    lp0 = n_mutation * log(p_mutation) + \
          (n_edges - n_mutation) * log(1 - p_mutation)
    return lp0


def _max_likelihood_tree_worker(tree_g, mc_df, cov_df, max_mutation=2, p_mutation=0.1, sub_tree_cutoff=12):
    top_n = 1

    n_edges = len(tree_g.edges)
    max_mutation = min(n_edges, max_mutation)

    record_names = mc_df.index

    if n_edges > sub_tree_cutoff:
        # cut the tree into left and right in the edge that has biggest betweenness_centrality
        # calculate best patterns for left and right separately, and then joint consider the overall pattern
        left_tree, right_tree, removed_edge = cut_by_highest_betweenness_centrality(tree_g)
        left_best_patterns, _ = _max_likelihood_tree_worker(
            left_tree,
            mc_df=mc_df.loc[:, mc_df.columns.isin(left_tree.nodes)],
            cov_df=cov_df.loc[:, cov_df.columns.isin(left_tree.nodes)],
            max_mutation=max_mutation, p_mutation=p_mutation, sub_tree_cutoff=sub_tree_cutoff)
        right_best_patterns, _ = _max_likelihood_tree_worker(
            right_tree,
            mc_df=mc_df.loc[:, mc_df.columns.isin(right_tree.nodes)],
            cov_df=cov_df.loc[:, cov_df.columns.isin(right_tree.nodes)],
            max_mutation=max_mutation, p_mutation=p_mutation, sub_tree_cutoff=sub_tree_cutoff)

        # for each DMR, go through all possible combination of best left and right pattern,
        # when not exceed max_mutation, also consider whether should we add the removed edge or not
        best_pattern_final = {}
        likelihood_final = {}
        for record_name in record_names:
            _this_mc_df = mc_df.loc[[record_name]]
            _this_cov_df = cov_df.loc[[record_name]]

            left_patterns = list(left_best_patterns[record_name]) + [()]  # add empty choice
            right_patterns = list(right_best_patterns[record_name]) + [()]  # add empty choice
            middle_patterns = [[removed_edge], []]

            # list all possible combined patterns
            pattern_dict = {}
            for left_i, left_pattern in enumerate(left_patterns):
                for right_i, right_pattern in enumerate(right_patterns):
                    for middle_pattern in middle_patterns:
                        joint_pattern = (list(left_pattern) if len(left_pattern) != 0 else []) + (
                            list(right_pattern) if len(right_pattern) != 0 else []) + (
                                            list(middle_pattern) if len(middle_pattern) != 0 else [])
                        _n_mutation = len(joint_pattern)
                        if _n_mutation > max_mutation:
                            continue

                        _this_group_mc_df, _this_group_un_mc_df = parse_one_pattern(
                            tree_g, joint_pattern, _this_mc_df, _this_cov_df)

                        # calculate tree likelihood on current pattern for all DMR
                        dmr_tree_likelihood = log_proba_beta_binomial(
                            _this_mc_df, _this_cov_df, _this_group_mc_df, _this_group_un_mc_df).values.sum()
                        # add mutation prior to tree likelihood, save to records
                        lp0 = mutation_likelihood(_n_mutation, p_mutation, n_edges)
                        try:
                            pattern_dict[_n_mutation][tuple(joint_pattern)] = dmr_tree_likelihood + lp0
                        except KeyError:
                            pattern_dict[_n_mutation] = {tuple(joint_pattern): dmr_tree_likelihood + lp0}
            _this_final_pattern = []
            _this_final_likelihood = []
            for _n_mutation, _n_mutation_patterns in pattern_dict.items():
                if _n_mutation != 0:
                    _s = pd.Series(_n_mutation_patterns).sort_values(ascending=False)[:top_n]
                    _this_final_pattern += _s.index.tolist()
                    _this_final_likelihood += _s.tolist()
                else:
                    # empty pattern
                    _this_final_pattern += [()]
                    _this_final_likelihood += list(_n_mutation_patterns.values())

            best_pattern_final[record_name] = np.array(_this_final_pattern)
            likelihood_final[record_name] = np.array(_this_final_likelihood)
        return pd.Series(best_pattern_final), pd.Series(likelihood_final)

    else:
        records = {}
        mutation_patterns = {}
        for n_mutation in range(1, max_mutation + 1):
            # Prior probability of the mutations, which is same for each n_mutation
            lp0 = n_mutation * log(p_mutation) + \
                  (n_edges - n_mutation) * log(1 - p_mutation)

            # each epimutation is removing one edge from the graph
            # for N epimutation, the result graph contain N + 1 disconnected sub-graph
            for i, edges in enumerate(combinations(tree_g.edges, n_mutation)):
                # get a and b for beta-binomial dist
                group_mc_df, group_un_mc_df = parse_one_pattern(tree_g, edges, mc_df, cov_df)

                # calculate tree likelihood on current pattern for all DMR
                dmr_tree_likelihood = log_proba_beta_binomial(mc_df, cov_df,
                                                              group_mc_df, group_un_mc_df).sum(axis=1)
                # add mutation prior to tree likelihood, save to records
                records[(n_mutation, i)] = dmr_tree_likelihood + lp0
                mutation_patterns[(n_mutation, i)] = edges

        # records_df: each row is a DMR record, each column is a (n_mutation, mutation_pattern_idx)
        records_df = pd.DataFrame(records)
        # mutation_pattern_series, index is (n_mutation, mutation_pattern_idx), value is the actual mutation pattern
        mutation_pattern_series = pd.Series(mutation_patterns)

        def __get_row_best_patterns(_row):
            _row_best_patterns = []
            _row_best_likelihoods = []
            for group, sub_row in _row.groupby(_row.index.get_level_values(0)):
                # index is pattern id, value is likelihood
                selected_pattern = sub_row.sort_values(ascending=False)[:top_n]
                _row_best_patterns.append(mutation_pattern_series.loc[selected_pattern.index])
                _row_best_likelihoods.append(selected_pattern)
            return pd.concat(_row_best_patterns).values, pd.concat(_row_best_likelihoods).values

        # top_n candidate pattern for each n_mutation
        patten_dict = {}
        likelihood_dict = {}
        for record_name, row in records_df.iterrows():
            _best_patterns, _likelihoods = __get_row_best_patterns(row)
            patten_dict[record_name] = _best_patterns
            likelihood_dict[record_name] = _likelihoods
        return pd.Series(patten_dict), pd.Series(likelihood_dict)


class DMRLineage:
    def __init__(self,
                 linkage,
                 linkage_anno_list,
                 mc_df: pd.DataFrame,
                 cov_df: pd.DataFrame,
                 max_mutation=5,
                 p_mutation=0.1,
                 sub_tree_cutoff=12):
        # validate names
        mc_names = mc_df.columns
        cov_names = cov_df.columns
        anno_names = pd.Index(linkage_anno_list)
        union_names = mc_names & cov_names & anno_names
        size_set = {mc_names.size, cov_names.size, anno_names.size, union_names.size}
        try:
            assert len(size_set) == 1
        except AssertionError:
            raise ValueError('Names in mc_df, cov_df and linkage_anno_list is inconsistent.')

        mc_df = mc_df[linkage_anno_list].copy()
        cov_df = cov_df[linkage_anno_list].copy()

        self.linkage = linkage
        self.linkage_anno_list = linkage_anno_list
        self.dendrogram = dendrogram(linkage, no_plot=True, labels=linkage_anno_list)
        self.tree_g = linkage_to_graph(linkage)
        self.n_leaves = self.linkage.shape[0] + 1

        # annotations
        self.cluster_id_to_name = {i: name for i, name in enumerate(linkage_anno_list)}
        self.cluster_name_to_id = {name: i for i, name in enumerate(linkage_anno_list)}
        self.non_singleton_node_dict = extract_all_nodes(linkage, labels=linkage_anno_list)

        self.mc_df_labeled = mc_df
        self.cov_df_labeled = cov_df

        self.mc_df = self.mc_df_labeled.copy()
        self.mc_df.columns = self.mc_df.columns.map(self.cluster_name_to_id)
        self.cov_df = self.cov_df_labeled.copy()
        self.cov_df.columns = self.cov_df.columns.map(self.cluster_name_to_id)

        self.max_mutation = max_mutation
        self.p_mutation = p_mutation
        self.sub_tree_cutoff = sub_tree_cutoff

        self.best_choice = None
        self.best_likelihood = None
        return

    def fit(self, cpu=1):
        tree_g = self.tree_g

        chunk_size = int(min(self.mc_df.shape[0], 10))
        futures = {}
        with ProcessPoolExecutor(cpu) as executor:
            for chunk_start in range(0, self.mc_df.shape[0], chunk_size):
                _chunk_mc_df = self.mc_df.iloc[chunk_start:chunk_start + chunk_size, :].copy()
                _chunk_cov_df = self.cov_df.iloc[chunk_start:chunk_start + chunk_size, :].copy()
                future = executor.submit(_max_likelihood_tree_worker,
                                         tree_g=tree_g,
                                         mc_df=_chunk_mc_df,
                                         cov_df=_chunk_cov_df,
                                         max_mutation=self.max_mutation,
                                         p_mutation=self.p_mutation,
                                         sub_tree_cutoff=self.sub_tree_cutoff)
                futures[future] = chunk_start

        results_dict = {}
        for future in as_completed(futures):
            chunk_start = futures[future]
            result = future.result()
            results_dict[chunk_start] = result

        dmr_best_mutations_list = []
        dmr_likelihoods_list = []
        for chunk_start in range(0, self.mc_df.shape[0], chunk_size):
            _patterns, _likelihoods = results_dict[chunk_start]
            dmr_best_mutations_list.append(_patterns)
            dmr_likelihoods_list.append(_likelihoods)
        total_records = pd.DataFrame({'mutation': pd.concat(dmr_best_mutations_list),
                                      'likelihoods': pd.concat(dmr_likelihoods_list)})
        self.best_choice = total_records.apply(lambda row: row['mutation'][row['likelihoods'].argmax()], axis=1)
        self.best_likelihood = total_records['likelihoods'].apply(lambda i: i.max())
        return

    def get_dmr_parsimony(self, dmr_id):
        raw_rate = (self.mc_df / self.cov_df).loc[dmr_id]
        tree_g = self.tree_g
        epi_mutation = self.best_choice[dmr_id]

        _this_mc_df = self.mc_df.loc[[dmr_id]].copy()
        _this_cov_df = self.cov_df.loc[[dmr_id]].copy()

        a, b = parse_one_pattern(tree_g, epi_mutation, _this_mc_df, _this_cov_df)
        tree_rate = (a / (a + b)).loc[dmr_id]
        data = pd.DataFrame([raw_rate, tree_rate, self.cov_df.loc[dmr_id]],
                            index=['raw', 'tree', 'cov']).T
        data.index.name = 'cluster'
        data = data.reindex(self.dendrogram['leaves'])
        data.index = self.dendrogram['ivl']
        data.reset_index(inplace=True)
        return data, epi_mutation

    def dmr_mutation_profile(self, dmr_id):
        _this_mc_df = self.mc_df.loc[[dmr_id]].copy()
        _this_cov_df = self.cov_df.loc[[dmr_id]].copy()

        _this_mc_df_with_name = self.mc_df_labeled.loc[[dmr_id]].copy()
        _this_cov_df_with_name = self.cov_df_labeled.loc[[dmr_id]].copy()
        epi_mutation = self.best_choice[dmr_id]

        mutation_profile = {}
        for edge in epi_mutation:
            parent = max(edge)
            this_child = min(edge)

            # find brother
            left, right, *_ = self.linkage.iloc[parent - self.n_leaves]
            left = int(left)
            right = int(right)

            _this_nodes_dict = self.non_singleton_node_dict[parent]
            _left_nodes = _this_nodes_dict['left']
            _right_nodes = _this_nodes_dict['right']

            if left == this_child:
                this_brother = right
                this_nodes = _left_nodes
                brother_nodes = _right_nodes
            elif right == this_child:
                this_brother = left
                this_nodes = _right_nodes
                brother_nodes = _left_nodes
            else:
                raise ValueError('Child node not in the parent linkage records')

            this_rate = _this_mc_df_with_name[this_nodes].values.sum() / _this_cov_df_with_name[this_nodes].values.sum()
            brother_rate = _this_mc_df_with_name[brother_nodes].values.sum() / _this_cov_df_with_name[
                brother_nodes].values.sum()

            mutation_profile[this_child] = this_rate - brother_rate
            mutation_profile[this_brother] = brother_rate - this_rate
        return pd.Series(mutation_profile, name=dmr_id)

    def mutation_delta(self, cluster):
        this_tree = self.tree_g.copy()
        data, mutations = self.get_dmr_parsimony(cluster)
        this_tree.remove_edges_from(mutations)

        nodes_group_rate_dict = {}
        for sub_g in nx.connected_component_subgraphs(this_tree):
            members = [i for i in sub_g.nodes if i < self.n_leaves]
            total_mc = self.mc_df.loc[cluster, members].sum()
            total_cov = self.cov_df.loc[cluster, members].sum()
            group_rate = total_mc / total_cov
            for member in sub_g.nodes:
                nodes_group_rate_dict[member] = group_rate

        results = {}
        for mutation_pair in mutations:
            child, parent = sorted(mutation_pair)
            mc_delta = nodes_group_rate_dict[child] - nodes_group_rate_dict[parent]
            results[child] = mc_delta
        return results, nodes_group_rate_dict

    def plot_dmr_lineage(self, dmr_id,
                         node_sizes=(10, 120),
                         node_size_norm=(0.1, 0.3),
                         hue_norm=(0.3, 0.8),
                         palette='viridis',
                         labelsize=7,
                         tree_bar_cax_tuple=None,
                         plot_node_id=False):

        if tree_bar_cax_tuple is None:
            fig = plt.figure(figsize=(6, 4), dpi=300)
            gs = fig.add_gridspec(20, 20, hspace=2, wspace=1)
            ax_tree = fig.add_subplot(gs[:15, :-1])
            ax_bar = fig.add_subplot(gs[15:, :-1])
            cax = fig.add_subplot(gs[15:, -1:])
        else:
            ax_tree, ax_bar, cax = tree_bar_cax_tuple

        # prepare plot data
        data, mutation = self.get_dmr_parsimony(dmr_id)
        mutation_profile, total_group_rate = self.mutation_delta(dmr_id)

        node_pos = plot_dendrogram(
            self.linkage,
            self.mc_df.columns,
            ax=ax_tree,
            plot_node_id=plot_node_id,
            node_palette=palette,
            node_hue=total_group_rate,
            node_hue_norm=hue_norm,
            node_size={k: abs(v)
                       for k, v in mutation_profile.items()},
            node_size_norm=node_size_norm,
            line_hue=total_group_rate,
            line_hue_norm=hue_norm,
            sizes=node_sizes)

        ymax = self.linkage[2].max()
        for k, v in mutation_profile.items():
            ax_tree.text(node_pos[k][0],
                         node_pos[k][1] + ymax * 0.05,
                         f'{v:.2f}',
                         fontsize=labelsize,
                         ha='left',
                         va='center',
                         rotation=90,
                         rotation_mode='anchor',
                         bbox=dict(
                             boxstyle="round",
                             linewidth=0.5,
                             fc=(1, 1, 1, 0.9),
                             ec=(0, 0, 0, 0.1),
                         ))

        ax_tree.set(xticks=[])
        sns.despine(ax=ax_tree, bottom=True)
        sns.despine(ax=ax_bar)

        plot_parsimony_data(data=data, ax=ax_bar, hue_norm=hue_norm, palette=palette)
        ax_bar.xaxis.set_tick_params(rotation=90)
        ax_bar.set(ylim=(0, 1))
        ax_bar.grid(axis='x', linewidth=0.3, alpha=0.5)
        ax_bar.set_xlim(ax_tree.get_xlim())

        plot_colorbar(cax,
                      cmap=palette,
                      cnorm=Normalize(*hue_norm),
                      hue_norm=hue_norm,
                      linewidth=0.7)

        ax_tree.set(ylabel='Cluster Distance')
        ax_bar.set(ylabel='mCG%', xlabel='Clusters')
        return ax_tree, ax_bar, cax

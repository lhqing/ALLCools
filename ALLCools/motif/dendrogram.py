from itertools import combinations

import networkx as nx
import numpy as np
from numpy import log
from scipy.special import betaln


def extract_all_nodes(linkage, labels=None):
    cluster_dict = {}
    cur_cluster_id = len(linkage) + 1
    total_obs = len(linkage) + 1
    if labels is None:
        labels = list(range(total_obs))
    for left, right, distance, n_obs in linkage:
        left = int(left)
        right = int(right)
        n_obs = int(n_obs)

        cluster_dict[cur_cluster_id] = {'left': set(),
                                        'right': set()}
        if n_obs == 2:
            left = labels[left]
            right = labels[right]
            # merge of 2 original observations
            cluster_dict[cur_cluster_id]['left'].add(left)
            cluster_dict[cur_cluster_id]['right'].add(right)
        else:
            # left and/or right are cluster
            if left < total_obs:
                left = labels[left]
                cluster_dict[cur_cluster_id]['left'].add(left)
            else:
                # node are cluster
                cluster_dict[cur_cluster_id]['left'].update(cluster_dict[left]['left'])
                cluster_dict[cur_cluster_id]['left'].update(cluster_dict[left]['right'])
            if right < total_obs:
                right = labels[right]
                cluster_dict[cur_cluster_id]['right'].add(right)
            else:
                # node are cluster
                cluster_dict[cur_cluster_id]['right'].update(cluster_dict[right]['left'])
                cluster_dict[cur_cluster_id]['right'].update(cluster_dict[right]['right'])
        cur_cluster_id += 1
    return cluster_dict


def dmr_distance(a: np.ndarray, b: np.ndarray):
    overlap = (a == b)[(a * b) != 0]
    return 1 - (overlap.sum() / overlap.size)


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


def log_proba_beta_binomial(x, n, a, b):
    """log likelihood for the beta-binomial dist, ignore part not related to a and b."""
    like = betaln((a + x), (b + n - x)) - betaln(a, b)
    return like


def parse_one_pattern(tree_g, edges_to_remove, mc_df, cov_df):
    """
    for a paticular epimutation combination (edges_to_remove),
    calculate the a and b for beta-binomial dist in each leaf node group.

    after removing the edges (epimutations),
    the leaf node group are leaf nodes in each of the disconnected sub graph.
    """
    group_mc_df = mc_df.copy()
    group_un_mc_df = cov_df - group_mc_df

    sub_g = tree_g.copy()
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


def max_likelihood_tree(linkage, mc_df, cov_df, max_mutation=2, p_mutation=0.1):
    tree_g = linkage_to_graph(linkage)
    n_edges = len(tree_g.edges)

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
    # select max tree for each item
    return records, mutation_patterns

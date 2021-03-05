import numpy as np
from sklearn.decomposition import PCA
from ..mcds.utilities import calculate_posterior_mc_frac
from scipy.stats import ks_2samp
from pynndescent import NNDescent
from sklearn.metrics import roc_curve
from random import choices
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns


class MethylScrublet:
    def __init__(self,
                 sim_doublet_ratio=2.0,
                 n_neighbors=None,
                 expected_doublet_rate=0.1,
                 stdev_doublet_rate=0.02,
                 metric='euclidean',
                 random_state=0,
                 n_jobs=-1):
        # initialize matrices
        self._mc_obs = None
        self._cov_obs = None
        self.n_obs = None
        self._frac_obs = None
        self._frac_sim = None
        self._pcs_obs = None
        self._pcs_sim = None
        self.clusters = None

        # initialize doublets score
        self.doublet_scores_obs_ = None
        self.doublet_scores_sim_ = None
        self.doublet_errors_obs_ = None
        self.doublet_errors_sim_ = None

        # Scrublet parameters
        self.sim_doublet_ratio = sim_doublet_ratio
        self.n_neighbors = n_neighbors
        self.expected_doublet_rate = expected_doublet_rate
        self.stdev_doublet_rate = stdev_doublet_rate
        self.random_state = random_state
        self.metric = metric
        self.n_jobs = n_jobs

        # doublets results
        self.predicted_doublets_ = None
        self.z_scores_ = None
        self.threshold_ = None
        self.detected_doublet_rate_ = None
        self.detectable_doublet_fraction_ = None
        self.overall_doublet_rate_ = None

    # Core Scrublet functions
    def fit(self, mc, cov, clusters=None):
        if isinstance(mc, xr.DataArray) and isinstance(cov, xr.DataArray):
            self._xarray_input = True
        elif isinstance(mc, np.ndarray) and isinstance(cov, np.ndarray):
            self._xarray_input = False
        else:
            raise TypeError('mc and cov should be both xr.DataArray or np.ndarray')

        np.random.seed(self.random_state)
        # calculate input posterior mC rate
        self._mc_obs = mc
        self._cov_obs = cov
        self.n_obs = mc.shape[0]
        self.clusters = clusters
        if self.n_neighbors is None:
            self.n_neighbors = min(50, int(round(0.5 * np.sqrt(self._mc_obs.shape[0]))))
        print('Calculating mC frac of observations...')
        self._frac_obs = calculate_posterior_mc_frac(mc_da=self._mc_obs,
                                                     cov_da=self._cov_obs)

        print('Simulating doublets...')
        self.simulate_doublets()

        print('PCA...')
        self.pca()

        print('Calculating doublet scores...')
        self.calculate_doublet_scores()
        self.call_doublets()
        return self.doublet_scores_obs_, self.predicted_doublets_

    def simulate_doublets(self):
        """ Simulate doublets by adding the counts of random observed transcriptome pairs."""
        n_sim = int(self.n_obs * self.sim_doublet_ratio)
        if self.clusters is not None:
            print('Cell cluster labels are given, will sample similar number of cells from each cluster.')
            clusters = self.clusters.reset_index(drop=True)
            major = clusters.value_counts()[0]
            use_cells = []
            for _, sub_series in clusters.groupby(clusters):
                replace = sub_series.size < major
                use_cells += pd.DataFrame(sub_series).sample(major, replace=replace).index.tolist()
            pair_ix = np.array([choices(use_cells, k=n_sim),
                                choices(use_cells, k=n_sim)]).T
        else:
            pair_ix = np.random.randint(0, self.n_obs, size=(n_sim, 2))

        # calculate mc frac of simulated data
        mc1 = self._mc_obs[pair_ix[:, 0], :]
        mc2 = self._mc_obs[pair_ix[:, 1], :]
        cov1 = self._cov_obs[pair_ix[:, 0], :]
        cov2 = self._cov_obs[pair_ix[:, 1], :]
        if self._xarray_input:
            mc1.coords['cell'] = range(n_sim)
            mc2.coords['cell'] = range(n_sim)
            cov1.coords['cell'] = range(n_sim)
            cov2.coords['cell'] = range(n_sim)
        self._frac_sim = calculate_posterior_mc_frac(mc1 + mc2, cov1 + cov2)
        return

    def pca(self):
        obs = self._frac_obs
        sim = self._frac_sim
        pca = PCA(n_components=min(200, obs.shape[1])).fit(obs)
        obs_pcs = pca.transform(obs)
        sim_pcs = pca.transform(sim)
        i = 0
        for i in range(obs_pcs.shape[1] - 1):
            cur_pc = obs_pcs[:, i]
            next_pc = obs_pcs[:, i + 1]
            p = ks_2samp(cur_pc, next_pc).pvalue
            if p > 0.01:
                break
        n_components = max(i + 1, 2)
        self._pcs_obs = obs_pcs[:, :n_components].copy()
        self._pcs_sim = sim_pcs[:, :n_components].copy()
        return

    def get_knn_graph(self, data):
        nn = NNDescent(data,
                       metric='euclidean',
                       n_jobs=self.n_jobs,
                       random_state=self.random_state)
        indices, distances = nn.query(data, k=self.n_neighbors + 1)
        knn = indices[:, 1:]
        return knn

    def calculate_doublet_scores(self):
        total_pcs = np.vstack((self._pcs_obs, self._pcs_sim))
        n_obs = self._pcs_obs.shape[0]
        n_sim = self._pcs_sim.shape[0]
        doublet_labels = np.concatenate((np.zeros(n_obs), np.ones(n_sim))).astype(int)

        # Adjust k (number of nearest neighbors) based on the ratio of simulated to observed cells
        k_adj = int(round(self.n_neighbors * (1 + n_sim / n_obs)))

        # Find k_adj nearest neighbors
        knn = self.get_knn_graph(total_pcs)

        # Calculate doublet score based on ratio of simulated cell neighbors vs. observed cell neighbors
        doublet_neighbor_mask = doublet_labels[knn] == 1  # get the identities of neighbors
        n_sim_neigh = doublet_neighbor_mask.sum(1)

        rho = self.expected_doublet_rate
        se_rho = self.stdev_doublet_rate
        r = self.sim_doublet_ratio
        nd = n_sim_neigh
        n = k_adj

        # Bayesian
        q = (nd + 1) / (n + 2)
        ld = q * rho / r / (1 - rho - q * (1 - rho - rho / r))
        se_q = np.sqrt(q * (1 - q) / (n + 3))
        se_ld = q * rho / r / (1 - rho - q * (1 - rho - rho / r)) ** 2 * np.sqrt(
            (se_q / q * (1 - rho)) ** 2 + (se_rho / rho * (1 - q)) ** 2)
        self.doublet_scores_obs_ = ld[doublet_labels == 0]
        self.doublet_scores_sim_ = ld[doublet_labels == 1]
        self.doublet_errors_obs_ = se_ld[doublet_labels == 0]
        self.doublet_errors_sim_ = se_ld[doublet_labels == 1]
        return

    def call_doublets(self, threshold=None):
        ld_obs = self.doublet_scores_obs_
        ld_sim = self.doublet_scores_sim_

        if threshold is None:
            y_score = np.concatenate([ld_obs, ld_sim])
            y_true = np.concatenate([np.ones_like(ld_obs), np.zeros_like(ld_sim)])
            fpr, tpr, thresholds = roc_curve(y_true, y_score)
            threshold = thresholds[np.argmin(tpr ** 2 + (1 - fpr) ** 2)]
            print(f'Automatically set threshold to {threshold:.2f}')

        se_obs = self.doublet_errors_obs_
        z = (ld_obs - threshold) / se_obs
        self.predicted_doublets_ = ld_obs > threshold
        self.z_scores_ = z
        self.threshold_ = threshold
        self.detected_doublet_rate_ = (ld_obs > threshold).sum() / len(ld_obs)  # FPR
        self.detectable_doublet_fraction_ = (ld_sim > threshold).sum() / len(ld_sim)  # TPR
        self.overall_doublet_rate_ = self.detected_doublet_rate_ / self.detectable_doublet_fraction_

        print(f'Detected doublet rate = {100 * self.detected_doublet_rate_:.1f}%')
        print(f'Estimated detectable doublet fraction = {100 * self.detectable_doublet_fraction_:.1f}%')
        print('Overall doublet rate:')
        print(f'\tExpected   = {100 * self.expected_doublet_rate:.1f}%')
        print(f'\tEstimated  = {100 * self.overall_doublet_rate_:.1f}%')
        return self.predicted_doublets_

    def plot(self):

        fig, (ax_roc, ax_hist) = plt.subplots(figsize=(6, 3), dpi=300, ncols=2)

        ax = ax_roc
        ld_obs = self.doublet_scores_obs_
        ld_sim = self.doublet_scores_sim_
        y_score = np.concatenate([ld_obs, ld_sim])
        y_true = np.concatenate([np.ones_like(ld_obs), np.zeros_like(ld_sim)])
        # tpr, fpr is reversed, because smaller score is better
        tpr, fpr, thresholds = roc_curve(y_true, y_score)
        sns.lineplot(ax=ax, x=fpr, y=tpr)
        ax.scatter([self.detected_doublet_rate_], [self.detectable_doublet_fraction_], color='red', zorder=99)
        ax.text(0.5, 0.1, f'Obs doub: {100 * self.detected_doublet_rate_:.1f}%\n'
                          f'Sim doub: {100 * self.detectable_doublet_fraction_:.1f}%',
                transform=ax.transAxes)
        ax.set(xlabel='TPR', ylabel='FPR', title='ROC')

        ax = ax_hist
        sns.histplot(self.doublet_scores_obs_, color='steelblue', ax=ax, label='Observed')
        sns.histplot(self.doublet_scores_sim_, color='gray', ax=ax, label='Simulated')
        ax.set(yticks=[], ylabel='', xlabel='Doublet Score')
        ax.legend()
        ymin, ymax = ax.get_ylim()
        ax.vlines(self.threshold_, ymin, ymax, color='red', linestyle='--')
        sns.despine(fig=fig)

        fig.tight_layout()
        return

    def _plot_cluster_dist(self):
        plot_data = pd.DataFrame({'groups': self.clusters,
                                  'Doublet Score': self.doublet_scores_obs_,
                                  'Is Doublet': self.predicted_doublets_})
        fig, ax = plt.subplots(figsize=(6, 2), dpi=300)
        sns.violinplot(data=plot_data, x='groups', y='Doublet Score', linewidth=0, scale='width')
        sns.stripplot(data=plot_data, x='groups', y='Doublet Score', s=1,
                      hue='Is Doublet', palette={True: 'red', False: 'black'})
        xmin, xmax = ax.get_xlim()
        ax.legend_.remove()
        ax.xaxis.set_tick_params(rotation=90)
        ax.hlines(self.threshold_, xmin, xmax, linestyle='--', linewidth=1, color='gray')
        sns.despine(ax=ax)

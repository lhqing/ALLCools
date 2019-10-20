import pandas as pd
from sklearn.cluster import MiniBatchKMeans


def dmr_kmeans(dmr_df, k, random_state=1):
    clf = MiniBatchKMeans(n_clusters=k,
                          init='k-means++',
                          max_iter=100,
                          batch_size=100,
                          verbose=0,
                          compute_labels=True,
                          random_state=random_state,
                          tol=0.0,
                          max_no_improvement=10,
                          init_size=None,
                          n_init=10,
                          reassignment_ratio=0.01)
    y_pred = clf.fit_predict(dmr_df.values)
    y_pred = pd.Series(y_pred, dmr_df.index)
    cluster_median = dmr_df.groupby(y_pred).median()
    return y_pred, cluster_median

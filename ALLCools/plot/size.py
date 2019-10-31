import seaborn as sns
from matplotlib.colors import Normalize

from .utilities import smart_number_format


def plot_sizebar(size_norm, sizes, ax):
    ticks = [0, 0.33, 0.66, 1]
    n_dots = len(ticks)

    snorm = Normalize(*size_norm, clip=True)

    def actual_size(i):
        return sizes[0] + snorm(i) * (sizes[1] - sizes[0])

    delta_norm = max(size_norm) - min(size_norm)
    actual_value = [min(size_norm) + i * delta_norm for i in ticks]
    actual_sizes = [actual_size(i) for i in actual_value]

    width = 1 / (n_dots + 1)

    x, y = zip(*[(0.5, width * (i + 1)) for i in range(n_dots)])
    ax.scatter(x, y, s=actual_sizes)
    sns.despine(left=True, bottom=True)
    ax.yaxis.tick_right()
    ax.set(xticks=[], yticks=y,
           yticklabels=[smart_number_format(i) for i in ticks])
    return ax

# Installation

## Setup Analysis Environment

### Install conda from miniconda or anaconda

To avoid potential conflicts with other packages, we strongly recommend you to use
a [conda environment](https://www.anaconda.com/products/individual). If you do not have a working installation or Python
3.6 (or later), consider installing [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (
see [Installing Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)).

Using such an isolated environment makes it possible to install a specific version of ALLCools with pip or conda and its
dependencies independently of any previously installed Python packages.

Note that you should always remember to activate the environment of your choice prior to running any Python command
whenever you start a new terminal session.

### Conda init

After installing conda, use conda init on your favorite shell.

```shell
# e.g., bash or zsh
conda init bash
# you need to restart shell after conda init
```

### Add channels

Run the following commands in their exact order to add bioconda into your conda channel.

```shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Create environment with required packages

First, you can create a `allcools_env.yaml` file as follows that contains the detail about the environment.

```yaml
name: allcools
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.8
  - pip
  - anndata
  - biopython
  - dask
  - numba
  - htslib>=1.9
  - jupyter
  - jupyter_contrib_nbextensions
  - leidenalg
  - natsort
  - netCDF4
  - networkx
  - opentsne
  - plotly
  - pybedtools
  - pyBigWig
  - pynndescent
  - pysam
  - pytables
  - scanpy
  - scikit-learn
  - seaborn
  - statsmodels
  - xarray
  - yaml
  - zarr
  - pip:
      - papermill
      - imblearn
      - allcools
```

Then, you can use the following command to create a conda environment called "allcools" and install all the required packages for you.

```shell
# first, install a package manager called mamba in base environment.
conda install mamba -n base -c conda-forge
# enter base env
conda activate base
# second, create a new environment and install all packages in the yaml file (content below)
mamba env create -f allcools_env.yaml
```

````{tip}
[mamba](https://mamba.readthedocs.io/en/latest/installation.html) is a CLI tool to manage conda environments. mamba can be installed alongside <code>conda</code> and it can provide faster sovles for big environments.  

We highly recommend you to use mamba for installing packages. If you don't want to use mamba, here is the conda command that achieves the same goal:
```shell
conda env create -f allcools_env.yaml
```
````

## Activate Environment

````{margin}
```{caution}
Remember to activate the "allcools" environment every time before you start any analysis.
```
````

```shell
# enter env
conda activate allcools

# exit env
conda deactivate
```

### Install optional packages

Here are some optional packages which might be hard to install on some old systems.

- <code>rpy2</code> (R and the R package pvclust) is used for the cluster dendrogram.
- <code>tpot</code> is used in REPTILE model.

```shell
mamba install -n allcools rpy2
mamba install -n allcools tpot xgboost dask-ml scikit-mdr skrebate
```

````{note}
If you cannot install rpy2, you can instead calculate a normal dendrogram with scipy function. In addition, if you cannot install tpot, REPTILE will automatically resort to the sklearn model. Note that <code>conda install</code> also works, but much slower than mamba. 
````

## Update

For updating the ALLCools package, you can enter the following commands by first activating your conda environment.

```shell
# enter env
conda activate allcools

# update package
pip install --upgrade allcools
```

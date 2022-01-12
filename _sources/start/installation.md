# Installation

## Setup Conda and Analysis Environment

### Install conda from miniconda or anaconda.

- IMPORTANT: select python 3
- [miniconda](https://docs.conda.io/en/latest/miniconda.html) (recommended)
- [anaconda](https://www.anaconda.com/products/individual) (larger)

### Conda init

After installed conda, use conda init on your favorite shell

```shell
# e.g., bash or zsh
conda init bash
# you need to restart shell after conda init
```

### Add channels

```shell
# run these command to add bioconda into your conda channel, 
# the order of these 3 line matters
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Create Environment With Required Packages

This command will create a conda environment called "allcools" and install all required packages for you. 
The `allcools_env.yaml` contains the detail about the environment, you can copy the content of that file below and create one by yourself.

```shell
# first, install a package manager called mamba in base environment
conda install mamba -n base -c conda-forge
# enter base env
conda activate base
# second, create a new environment and install all packages in the yaml file (content below)
mamba env create -f allcools_env.yaml
```

````{tip}
conda is very slow when solving a large number of packages, [mamba](https://mamba.readthedocs.io/en/latest/installation.html) is a lightening fast replacement of conda. I highly recommend you use mamba to install packages. If you don't want to use mamba, here is the conda command that achieves the same goal:
```shell
conda env create -f allcools_env.yaml
```
````

### Content of `allcools_env.yaml`

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
  - jupyter
  - jupyter_contrib_nbextensions
  - leidenalg
  - natsort
  - networkx
  - opentsne
  - plotly
  - pybedtools
  - pyBigWig
  - pynndescent
  - pysam
  - scanpy
  - scikit-learn
  - seaborn
  - statsmodels
  - xarray
  - yaml
  - pip:
    - papermill
    - imblearn
    - allcools
```

## Activate Environment
````{margin}
```{caution}
Remember to enter the environment every time before you start any analysis.
```
````

```shell
# enter env
conda activate allcools

# exit env
conda deactivate
```


### Install optional packages
Here are some additional packages which might be hard to install for some old systems and are optional.
```shell
# rpy2 (R and the R package pvclust) is used for the cluster dendrogram, 
# if you can't install, will just calculate a normal dendrogram with scipy function
mamba install -n allcools rpy2
# tpot is used in REPTILE model, but if not found, REPTILE will turn to use normal sklearn model
mamba install -n allcools tpot xgboost dask-ml scikit-mdr skrebate
# note that conda install also works, just change mamba to conda. But its usually slower.
```

## Update

```shell
# enter env first
conda activate allcools
pip install --upgrade allcools
```

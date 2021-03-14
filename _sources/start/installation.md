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
The `allcools_env.yaml` contains the detail about the environment, you can copy the content of that file below.
```shell
# conda is a bit slow, this step will take ~30 min in my server, 
# but you only need to do this once.
conda env create -f allcools_env.yaml
```

#### Content of `allcools_env.yaml`
```yaml
name: allcools
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.7
  - pip
  - jupyter
  - jupyter_contrib_nbextensions
  - pandas
  - seaborn
  - matplotlib
  - pip:
    - papermill
# TODO this list is not completed yet
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

## Install ALLCools
### Install
```shell
# enter env first
pip install allcools
```

### Update
```shell
# enter env first
pip install --upgrade allcools
```

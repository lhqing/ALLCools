# ALLCools

- A toolkit for ALLC format and bulk/single-cell methylation sequencing analysis. 
- It provides all necessary functions about generating, manipulating and transfering 2 file formats:
  - the ALLC format, a tab-sparated text format storing **base-level** raw methylation information for **1 sample**.
  - the MCDS format, a integrated HDF5 dataset storing **region-level** raw methylation information for **a group of samples/cells**.
- **This package is still under active development**

## ALLCools Function Map
![](/doc/file/ALLCools.svg)

## Get Start
### requirements
ALLCools requires python > 3.6

### installation
```shell
git clone https://github.com/lhqing/ALLCools.git
cd ALLCools
pip install .
```

### update
```shell
# cd /path/to/ALLCools
git pull
pip install .
```

### use the ALLCools
- CLI (recommand for non-python users)
```shell
allcools -h
allcools [FUNCTION_NAME] -h
```
- API (recommand for integrating ALLCools with other python package, such as pybedtools, pybigwig, deeptools.)

```python
from ALLCools.api import *
```

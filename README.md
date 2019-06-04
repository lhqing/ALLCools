# ALLCools

- A toolkit for ALLC format and bulk/single-cell methylation sequencing analysis. 
- It provides all necessary functions about generating, manipulating and transfering the ALLC format, a simple tab-sparated text format storing base-level raw methylation information.
- **This package is still under active development**

## ALLCools Function Map
![](/doc/file/ALLCools.svg)

## Get Start
```shell
# requirements
# ALLCools requires python > 3.6

# installation
git clone https://github.com/lhqing/ALLCools.git
cd ALLCools
pip install .

# update (necessary)
# cd /path/to/ALLCools
git pull
pip install .

# use the ALLCools
# CLI (recommand for non-python users)
allcools -h
allcools [FUNCTION_NAME] -h
```

```python
# API (recommand for integrate ALLCools with other python package, 
# such as pybedtools, pybigwig, deeptools.)
from ALLCools.api import *

```

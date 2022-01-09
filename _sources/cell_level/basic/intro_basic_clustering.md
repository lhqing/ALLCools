# Basic Clustering Walk-through

## Section table of contents

```{tableofcontents}
```

## Input Files

### For 100Kb bins clustering
The dataset we used for 100Kb clustering documentation comes from two brain regions of adult P56 male mice. 
The MOp stand for primary motor cortex, the HIP stand for hippocampus (CA) region. 
These two datasets belong to a larger whole brain methylome atlas dataset (add reference).

#### Download Input Files
- Cell metadata: ADD DOWNLOAD URL
- single-cell ALLC files: ADD DOWNLOAD URL
- MCDS files: ADD DOWNLOAD URL

### For 5Kb bins clustering
The dataset we used for 5Kb clustering documentation comes from pituitary of adult P56 male mice (add reference).

#### Download Input Files
- Cell metadata: ADD DOWNLOAD URL
- single-cell ALLC files: ADD DOWNLOAD URL
- MCAD files: ADD DOWNLOAD URL

## Prepare your own datasets
### Cell metadata
Using the sample dataset above as an example. The cell metadata file is a user-defined file that has no fixed format. 
It usually contains mapping metric and experimental design metadata for each cell. 
We use this file to perform some basic filtering to remove low quality cells.

### Single-cell ALLC files
See [`allcools bam-to-allc`](../../command_line/allcools_allc.ipynb) for how to generate single-cell ALLC files 
from filtered single-cell BAM files.

### MCDS files
See [`allcools mcds`](../../command_line/allcools_mcds.ipynb) for how to generate MCDS files 
from single-cell ALLC files

### MCAD files
See [`allcools mcad`](../../command_line/allcools_mcad.ipynb) for how to generate MCAD files 
from single-cell ALLC files

# Basic Clustering Walk-through

## Section table of contents

```{tableofcontents}
```

## Input Files

### For 100Kb bins clustering
The dataset we used for 100Kb clustering documentation comes from the hippocampus (HIP) of adult P56 male mice. These datasets belong to a larger mouse brain methylome atlas dataset {cite:p}`Liu2021`.

#### Download Input Files
- Cell metadata: ADD DOWNLOAD URL
- single-cell ALLC files: ADD DOWNLOAD URL
- MCDS files: ADD DOWNLOAD URL

### For 5Kb bins clustering
The dataset we used for 5Kb clustering documentation comes from human PBMC (ADD REFERENCE).

#### Download Input Files
- Cell metadata: ADD DOWNLOAD URL
- single-cell ALLC files: ADD DOWNLOAD URL
- MCDS files: ADD DOWNLOAD URL

## Prepare your own datasets
### Cell metadata
Using the sample dataset above as an example. The cell metadata file is a user-defined file that has no fixed format. 
It usually contains mapping metric and experimental design metadata for each cell. 
We use this file to perform some basic filtering to remove low quality cells.

### Single-cell ALLC files
See [`allcools bam-to-allc`](../../command_line/allcools_allc.ipynb) for how to generate single-cell ALLC files 
from filtered single-cell BAM files.

### MCDS files
See [`allcools generate-dataset`](../../command_line/allcools_dataset.ipynb) for how to generate MCDS files containing raw count matrix, and/or hypo-/hyper-methylation score matrix for any kinds of genomic features.


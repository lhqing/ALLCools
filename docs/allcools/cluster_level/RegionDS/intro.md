# Post-clustering Genomic Analysis

```{tableofcontents}
```

## Prepare Cluster Level Files

### Merge Pseudo-bulk ALLC Files

After clustering analysis, we can merge single-cell ALLC files into pseudo-bulk ALLC files using 
[`allcools merge`](../../command_line/allcools_merge.ipynb). The merged ALLC files have the same format as single-cell ALLC files, except for the fifth column (mc) and sixth column (cov) are sum from single cells.

### Extract Certain Methylation Context

Differential Methylated Region (DMR) analysis will be focus on the CpG context, 
to extract all the CpG sites from the ALLC file, you can use the 
[`allcools extract`](../../command_line/allcools_extract.ipynb). 
The resulting file is still in the ALLC format.

Alternatively, if CpH sites are not of your interest, you can extract CpG sites from all the single-cell ALLC files before doing the merge ALLC step, as this will greatly fasilitate the merging step (CpG is ~20 times less than CpH in mm10).

#### Merge CpG strand

Merge the paired CpG sites could increase the performance of DMR calling algorithms, 
since paired CpG sites should have the same methylation status. 
To achieve this, just add an `--strandness merge` flat in the `allcools extract` command.

### Generate BigWig File From ALLC File

To generate BigWig file from ALLC file, you can use the 
[`allcools bw`](../../command_line/allcools_bw.ipynb). There are two kinds of BigWig files can be generated, 
the `frac.bw` stores the methylation fraction (mc/cov) of each entry, 
the `cov.bw` stores the coverage (cov) of each entry.
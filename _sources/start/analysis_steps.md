# Analysis Steps

ALLCools separates the analysis into two parts:
1. Cell-level analysis, including
    - Clustering analysis
    - Differentially Methylated Gene (DMG) analysis
    - Data integration
    - Identify potential doublets
2. Cluster-level analysis, including
    - Differentially Methylated Region (DMR) analysis
    - DMR annotation
    - DMR motif analysis (upstream regulators of DMRs)
    - DMR - Gene correlation analysis (downstream targets of DMRs)

In general, the **cell-level analysis** is focused on individual cells' overall diversity 
using relatively large features (kilobase-level regions or whole gene body). 
One of the key results from **cell-level analysis** is the identification of **cell clusters**.
The **cluster-level analysis** will use these cell cluster labels to merge the cells into pseudo-bulk methylomes. 
The high-coverage pseudo-bulk methylomes allow us to analyze methylation diversity at higher genome resolution 
(hundred bps or even single base).

## Cell Level Analysis

### Sections
- [Basic walk-through of the clustering analysis](../cell_level/basic/intro_basic_clustering.md).
- [Step-by-step description of the clustering analysis](../cell_level/step_by_step/intro_step_by_step_clustering.md).
- [Identification of cluster Differentially Methylated Genes](../cell_level/dmg/intro_dmg.md).
- [Cell-level data integration](../cell_level/integration/intro_integration.md).
- [Potential cell doublets identification](../cell_level/doublets/intro_doublets.md).

### Input
1. Cell metadata contains quality metrics (Defined by user, used in basic cell filtering).
2. Cell-by-feature raw count matrix (MCDS files) generated from single-cell ALLC files via [`allcools mcds`](
  ../command_line/allcools_mcds.ipynb) for large genomic regions (e.g., 100Kb bin) and gene body counts.
3. (AND/OR) Cell-by-feature hypo-methylation score matrix (MCAD files) generated from single-cell ALLC files via 
  [`allcools mcad`](../command_line/allcools_mcad.ipynb)

### Two Clustering Strategies for different tissue type

#### Clustering using raw counts from 100Kb genomic bins
In this strategy, we start from the cell-by-100kb-bin raw count matrix ([MCDS](mcds-fig)) for clustering analysis. 
A quick demo is [here](../cell_level/basic/mch_mcg_100k_basic.ipynb), the step-by-step descriptions can be found 
[here](../cell_level/step_by_step/100kb/intro_100kb.md).

We found {cite}`Liu2020-ab` using 100Kb genomic bins provides good clustering ability and is also computationally 
efficient due to the relatively small number of features (27K features in mm10, v.s. 540K features using 5Kb bins). 
We also tested the effect of different feature sizes in clustering adult mouse brain single-cell methylomes, 
and found the size of features do not have major impact on the clustering results  {cite}`Liu2020-ab`. 
However, this conclusion may not apply to all tissues, as the scale of methylation diversity might be different in 
different tissue or cell types. We suggest using the 100Kb bin size as a start point, and test other bin sizes when 
necessary. The [`allcools mcds`](../command_line/allcools_mcds.ipynb) does provide the flexibility to 
generate feature matrix at any bin size or through user-defined region sets.

#### Clustering using hypo-methylation score from 5Kb genomic bins
We also noticed that for some non-brain tissues (pituitary, PBMC etc.), the 
100Kb bin clustering strategy had a hard time identifying some known cell types. One possible explanation is that the 
methylation diversity of these cell types mainly occur at small discontinuous regulatory regions (DMRs) while the large 
100Kb bins can only capture DMV or large-hypo DMR level diversities. 
To solve this problem, we developed another strategy by using 
[algorithms adapted from snATAC-seq](../cell_level/step_by_step/5kb/intro_5kb.md) 
on small genomic bins (5Kb by default).

Specifically, this strategy start from cell-by-5kb-bin hypo-methylation score matrix ([MCAD](mcad-fig)) for clustering 
analysis. A quick demo is [here](../cell_level/basic/mcg_5kb_basic.ipynb), the step-by-step descriptions can be found
[here](../cell_level/step_by_step/5kb/intro_5kb.md).

## Cluster Level Analysis

### Sections
- Prepare pseudo-bulk ALLC files
- Call Differentially Methylated Region (DMR)
- DMR annotation
- DMR motif analysis (finding upstream regulators of DMRs)
- DMR - Gene correlation analysis (finding downstream targets of DMRs)

### Pseudo-bulk Files
Starting from **single-cell ALLC files** and **cluster labels** from the cell-level analysis, there are three steps to 
prepare cluster-level data files.

1. We create pseudo-bulk ALLC files by merging from single-cell ALLC files via 
[`allcools merge`](../command_line/allcools_merge.ipynb). 

2. We then use proper algorithms ([see discussion](../discuss/dmr_benchmark.md)) to identify DMRs. 
The expected results of this step are a DMR region list and associated statistics.

3. Using the DMR regions and additional datasets, we can annotate DMRs with multiple kinds of information, and store the 
DMR-by-feature information in anndata.AnnData format, where the `obs` dimension are the DMRs, the `var` dimension are 
all kinds of features.

The data model is illustrated in the following figure:

```{figure} ./cluster-level-analysis.png
---
height: 300px
name: cluster-level-analysis-fig
---
Cluster level analysis data model.
```

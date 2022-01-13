# Analysis Steps

ALLCools covers two parts of analysis:

1. Cellular analysis, including
   - [Clustering analysis](../cell_level/basic/intro_basic_clustering)
   - [Differentially Methylated Gene (DMG) analysis](../cell_level/dmg/intro_dmg)
   - [Data integration](../cell_level/integration/intro_integration)
   - [Identify potential doublets](../cell_level/doublets/intro_doublets)
2. Genomic analysis, including
   - [Differentially Methylated Region (DMR) analysis](../cluster_level/RegionDS/intro)
   - [DMR annotation](../cluster_level/RegionDS/02.annotation)
   - [DMR motif analysis (upstream regulators of DMRs)](../cluster_level/RegionDS/intro_motif)
   - [DMR - Gene correlation analysis (downstream targets of DMRs)](../cluster_level/Correlation/intro_corr)
   - [Enhancer prediction](../cluster_level/REPTILE/intro_reptile)

In general, the **cellular analysis** is focused on individual cells' overall diversity using relatively large features (kilobase-level regions or whole gene body). One of the key results from **cellular analysis** is the identification of **cell clusters**. We will then merge the single-cell methylome and other profiles to pseudo-bulk level based on the cluster labels, which gives us enough coverage to perform **genomic analysis** at hundred-bp or even single-base resolution.

## Cellular Analysis

### Sections
- [Basic walk-through of the clustering analysis](../cell_level/basic/intro_basic_clustering.md).
- [Step-by-step description of the clustering analysis](../cell_level/step_by_step/intro_step_by_step_clustering.md).
- [Identification of  Differentially Methylated Genes clusters](../cell_level/dmg/intro_dmg.md).
- [Cell-level data integration](../cell_level/integration/intro_integration.md).
- [Potential cell doublets identification](../cell_level/doublets/intro_doublets.md).


### Input

1. Cell metadata contains quality metrics (defined by user; used in basic cell filtering).
2. Cell-by-feature raw count matrix (MCDS files) generated from single-cell ALLC files via [`allcools mcds`](
   ../command_line/allcools_mcds.ipynb) for large genomic regions (e.g., 100Kb bin) and gene body counts.
3. (AND/OR) Cell-by-feature hypo-methylation score matrix (MCAD files) generated from single-cell ALLC files via 
   [`allcools mcad`](../command_line/allcools_mcad.ipynb).

### Two Clustering Strategies for different tissue types

#### Clustering using raw counts from 100Kb genomic bins

In this strategy, we start from the cell-by-100kb-bin raw count matrix ([MCDS](mcds-fig)) for clustering analysis. 
 [Here](../cell_level/basic/mch_mcg_100k_basic.ipynb) is a quick demo and its step-by-step descriptions can be found 
[here](../cell_level/step_by_step/100kb/intro_100kb.md).

We found {cite}`Liu2021` using 100Kb genomic bins provides good clustering ability and is also computationally 
efficient due to the relatively small number of features (27K features in mm10, v.s. 540K features using 5Kb bins). 
We also tested the effect of different feature sizes when clustering adult mouse brain single-cell methylomes and found that feature sizes do not have a major impact on the clustering results  {cite}`Liu2021`. 
However, this conclusion may not apply to all tissues, as the scale of methylation diversity might be different in 
different tissues or cell types. We suggest using the 100Kb bin size as a starting point and testing other bin sizes when 
necessary. That being said, the [`allcools mcds`](../command_line/allcools_mcds.ipynb) does provide the flexibility to 
generate feature matrix at any bin size or through user-defined region sets.

#### Clustering using hypo-methylation score from 5Kb genomic bins

We also noticed that for some non-brain tissues (pituitary, PBMC etc.), the 
100Kb bin clustering strategy had a hard time identifying some known cell types. One possible explanation is that the 
methylation diversity of these cell types mainly occur at small discontinuous regulatory regions (DMRs) while the large 
100Kb bins can only capture DMG or large-hypo DMR level diversities. 
To solve this problem, we developed another strategy by using 
[algorithms adapted from snATAC-seq](../cell_level/step_by_step/5kb/intro_5kb.md) 
on small genomic bins (5Kb by default).

Specifically, this strategy starts from a cell-by-5kb-bin hypo-methylation score matrix ([MCAD](mcad-fig)) for clustering 
analysis. [Here](../cell_level/basic/mcg_5kb_basic.ipynb) is a quick demo and its step-by-step descriptions can be found
[here](../cell_level/step_by_step/5kb/intro_5kb.md).

## Genomic Analysis

### Sections

- [Prepare pseudo-bulk ALLC files](../cluster_level/RegionDS/intro)
- [Call Differentially Methylated Region (DMR)](../cluster_level/RegionDS/01a.call_dmr)
- [DMR annotation](../cluster_level/RegionDS/02.annotation.ipynb)
- [DMR motif analysis (finding upstream regulators of DMRs)](../cluster_level/RegionDS/intro_motif.md)
- [DMR - Gene correlation analysis (finding downstream targets of DMRs)](../cluster_level/Correlation/intro_corr)
- [Enhancer prediction with REPTILE algorithm](../cluster_level/REPTILE/intro_reptile.md)

### Pseudo-bulk Files

Starting from **single-cell ALLC files** and **cluster labels** from the cell-level analysis, there are three steps for preparing cluster-level data files.

1. We create pseudo-bulk ALLC files by merging from single-cell ALLC files via 
   [`allcools merge`](../command_line/allcools_merge.ipynb). 

2. We then use proper algorithms ([see discussion](../discuss/dmr_benchmark.md)) to identify DMRs. 
   The expected results of this step are a DMR region list and associated statistics.

3. Using the DMR regions and additional datasets, we can annotate DMRs with multiple kinds of information and store the 
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

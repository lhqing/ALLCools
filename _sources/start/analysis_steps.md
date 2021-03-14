# Analysis Steps

ALLCools separates the analysis into two parts:
1. Cell-level analysis, including
    - Clustering Analysis
    - Differentially Methylated Gene (DMG) Analysis
    - Data Integration
2. Cluster-level analysis, including
    - Differentially Methylated Region (DMR) Analysis
    - DMR Annotation
    - DMR Motif Analysis
    - DMR - Gene Correlation Analysis
    - Enhancer Prediction

In general, the **cell-level analysis** is focused on individual cells' overall diversity 
using relatively large features (kilobase-level regions or whole gene body). 
One of the key results from **cell-level analysis** is the identification of **cell clusters**.
The "cluster-level analysis" will use the cluster labels to merge pseudo-bulk methylomes from 
single-cell methylomes. The high-coverage pseudo-bulk methylomes allow us to analyze cluster-level 
methylation diversity at higher genome resolution (hundred bps or even single base).


## Cell Level Analysis
### Input Files


### Tow Clustering Strategies
#### Clustering using raw counts from 100Kb genomic bins


#### Clustering using hypo-methylation score from 5Kb genomic bins


## Cluster Level Analysis
### Input Files

### DMR Calling

### DMR Analysis

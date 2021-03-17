ALLCools: ALL Cytosine tools
============================

ALLCools is a python package for single-cell DNA methylome data analysis. 
It includes functions to prepare input files, perform cell-level analysis 
(e.g., clustering), and cluster-level analysis 
(e.g., DMR calling and following analysis).


## Documentation Organization
1. To install and get an overview, please read the **GET STARTED** part;
2. To generate analysis files, please read the **COMMAND LINE TOOLS** part;
3. The analysis steps in ALLCools are organized in two parts:
    1. The **CELL LEVEL ANALYSIS** part go through the clustering, cluster DMG, and integration analysis;
    2. The **CLUSTER LEVEL ANALYSIS** part go through the DMR and enhancer calling analysis;
4. To customize visualization, please read the **VISUALIZATION** part;
5. The **DISCUSSION** part provides more details on the analysis steps;
6. The **API** part provides reference to each function.

```{figure} ./doc_organize.png
---
height: 250px
name: doc-organize-fig
---
ALLCools documentation organization.
```


## Author
Hanqing Liu

## Support
- The source code is on [github](https://github.com/lhqing/ALLCools);
- For releases and changelog, see [github releases page](https://github.com/lhqing/ALLCools/releases);
- For bugs and feature requests, please use the [issue tracker](https://github.com/lhqing/ALLCools/issues).

## Citation
When using ALLCools, please cite our paper and related methods papers. 
For more details, please read the [citation page](./project_info/citation.md).

[Hanqing Liu, Jingtian Zhou, Wei Tian, Chongyuan Luo, Anna Bartlett, Andrew Aldridge, Jacinta Lucero, 
Julia K. Osteen, Joseph R. Nery, Huaming Chen, Angeline Rivkin, Rosa G Castanon, Ben Clock, Yang Eric Li, 
Xiaomeng Hou, Olivier B. Poirion, Sebastian Preissl, Antonio Pinto-Duarte, Carolyn O’Connor, Lara Boggeman, 
Conor Fitzpatrick, Michael Nunn, Eran A. Mukamel, Zhuzhu Zhang, Edward M. Callaway, Bing Ren, Jesse R. Dixon, 
M. Margarita Behrens, Joseph R. Ecker. 2020. 
“DNA Methylation Atlas of the Mouse Brain at Single-Cell Resolution.” bioRxiv.
](https://doi.org/10.1101/2020.04.30.069377)
ALLCools: ALL methyl-Cytosine tools
============================

With the advancement of genome- and methylome-sequencing, ALLCools is a software package for single-cell DNA methylome data analysis. This user-friendly tool supports a broad range of next-generation sequencing (NGS) studies, from cellular analysis to genomic analysis and more. 

## Documentation Organization

1. For installation, please read **GET STARTED** [here](https://lhqing.github.io/ALLCools/start/installation.html);
2. The ALLCools workflow encompasses the following two secitons:
   a. **CELLULAR ANALYSIS**: This [section](https://lhqing.github.io/ALLCools/cell_level/basic/intro_basic_clustering.html) goes over basic clustering steps, differential methylated genes (DMG) analsysis, data integration, and doublets identification;
   b. **GENOMIC ANALYSIS**: This [section](https://lhqing.github.io/ALLCools/cluster_level/intro.html) covers differential methylated region (DMR) calling, genome annotation, DNA motif analysis, correlation analysis, and enhancer prediction.
3. To generate analysis files from your own data, please read the **COMMAND LINE TOOLS** [here](https://lhqing.github.io/ALLCools/command_line/allcools.html);
4. To customize visualization, please read the **VISUALIZATION** [here](https://lhqing.github.io/ALLCools/visualization/intro.html);
5. For more detials on the analysis steps, please read the **DISCUSSION** [here](https://lhqing.github.io/ALLCools/discuss/intro.html);
6. For references to functions used in this package, please read the **API** [here](https://lhqing.github.io/ALLCools/api/ALLCools/index.html).

```{figure} ./img/doc_organize.png
---
height: 350px
name: doc-organize-fig
---
ALLCools documentation organization.
```
<span style="color:red">Note: this figure can be a bit confusing wrt how you position the command line tools and visualization chunks. 1) Are the files prepared in the coomand line tools also applicable to the genome analysis? 2) What does browser mean in visualization? </span>.


## Authors

- Hanqing Liu, developer, initial conception
- Jingtian Zhou, developer, 5kb clustering algorithms

## Support

```{figure} ./img/open_issue.png
---
height: 150px
name: open-issue-fig
figclass: margin
---
Click on this to create a page specific issue.
```
<span style="color:red">Note: this figure is wrongly positioned on the page. I think it's better to just get rid of it. Users will know where to find it. </span>.


- The source code is on [github](https://github.com/lhqing/ALLCools);
- For releases and changelog, please check out the [github releases page](https://github.com/lhqing/ALLCools/releases);
- For bugs and feature requests, please use the [issue tracker](https://github.com/lhqing/ALLCools/issues).
- For page-specific issues, please use the "open issue" button on the top-right toggle.

## Citing ALLCools

If ALLCools has been significant in your research, and you would like to acknowledge the project in your academic publication, we suggest citing our paper {cite:p}`Liu2021`. For any specific methods and algorithm, also consider citing the original author's paper included in the Citation and Reference page [here](https://lhqing.github.io/ALLCools/project_info/citation.html). 

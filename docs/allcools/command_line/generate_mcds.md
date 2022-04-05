# allcools generate-dataset

## Tutorial 

**Example 1:**

This example covers comprehensively the situation where five features denoted after regions (e.g. chrom100k, chrom5k, geneslop2k, promoter, CGI) and their corresponding contexts specified with count (e.g. CGN, CHN, CAN, CCC) are used for generating mcds files, while two features denoted after quantifiers (e.g. promoter, CGI) and their respective contexts specified with hypo/hyper-score (e.g. CGN) are used for generating mcad files within the dataset. 

    allcools generate-dataset 
    --allc_table allc_table.tsv
    --output_path Farlik2016CSC.mcds
    --chrom_size_path hg38.main.nochrM.chrom.sizes
    --obs_dim cell 
    --cpu 20
    --chunk_size 50
    --regions chrom100k 100000
    --regions chrom5k 5000
    --regions geneslop2k hg38_genecode_v28.geneslop2k.bed.gz
    --regions promoter promoter_slop2k.sorted.bed.gz
    --regions CGI hg38_CGI.sorted.bed.gz
    --quantifiers chrom100k count CGN, CHN, CAN, CCC
    --quantifiers chrom5k count CGN
    --quantifiers geneslop2k count CGN
    --quantifiers promoter count CGN
    --quantifiers CGI count CGN
    --quantifiers promoter hypo-score CGN cutoff=0.9
    --quantifiers CGI hypo-score CGN cutoff=0.9

**Notes:**

    --allc_table allc_table.tsv
Specify the absolute file path of the allc table file in this line. Here is an example of what the allc table looks like:

> GSM1922083_A01  /gale/netapp/cemba3c/projects/SingleCellMethylomes/Angermueller2016NM/allc/GSM1922083_A01.allc.tsv.gz

The first part indicates the cell name (e.g. GSM1922083_A01) whereas the second part indicates the allc file path of the cell (e.g. /gale/netapp/cemba3c/projects/SingleCellMethylomes/Angermueller2016NM/allc/GSM1922083_A01.allc.tsv.gz). Make sure the two parts are separated by a tab. 

    --output_path Farlik2016CSC.mcds
 Specify here the absolute path of the output file and how you'd like to name the file.  

    --obs_dim cell 
Use cell as the obs_dim. 

    --cpu 20
Specify here how much cpu usage you need for running the command. 

    --chunk_size 50
Determine the parallel chunk size. The default size is 50. 

    --regions chrom100k 100000

    --regions chrom5k 5000

    --regions geneslop2k hg38_genecode_v28.geneslop2k.bed.gz

    --regions promoter promoter_slop2k.sorted.bed.gz

    --regions CGI hg38_CGI.sorted.bed.gz
Specify the features and datasets needed for generating the MCDS files. Features include chrom100k, chrom5k, geneslop2k, promoter, CGI, etc. For the chrom100k and chrom5k, specify the bin size using 100000 and 5000 respectively. For geneslop2k, promoter and CGI, specify the absolute bed file paths. 

    --quantifiers chrom100k count CGN, CHN, CAN, CCC
    --quantifiers chrom5k count CGN
    --quantifiers geneslop2k count CGN
    --quantifiers promoter count CGN
    --quantifiers CGI count CGN
Specify here the mc context needed for each feature for generating the MCDS files. Use "count" quantifier to represent the data type of MCDS count matrix.

    --quantifiers promoter hypo-score CGN cutoff=0.9
Specify here he mc context needed for each feature for generating the MCAD files. Use hypo-score to indicate the probability of the bins that are hypo-methylated. Inversely, use hyper-score to indicate the probability of the bins that are hyper-methylated.  

**Example 2:**
This example shows a simpler case where three features denoted after regions (e.g. chrom100k, chrom5k, geneslop2k) and their corresponding contexts specified with count (e.g. CGN, CHN) are used for generating mcds files, while one features denoted after quantifiers (e.g. chrom5k) and its contexts specified with hypo-score (e.g. CGN) are used for generating mcad files within the dataset. 

    allcools generate-dataset   
    --allc_table Luo2017Science_Human_allc_table.tsv   
    --output_path Luo2017Science_Human.mcds   
    --chrom_size_path hg19.main.nochrM.chrom.sizes   
    --obs_dim cell   
    --cpu 20   
    --chunk_size 50   
    --regions chrom100k 100000   
    --regions chrom5k 5000   
    --regions geneslop2k gencode.v28lift37.annotation.gene.bed.gz
    --quantifiers chrom100k count CGN,CHN  
    --quantifiers geneslop2k count CGN,CHN  
    --quantifiers chrom5k hypo-score CGN cutoff=0.9  

 

**Example 3:**
This example is similar to the ones above. The major difference is that it includes a mixture of hypo and hyper-score to denote the value and reverse value of mc contexts. 

    allcools generate-dataset   
    --allc_table allc_table.tsv   
    --output_path Li2018NCB.mcds   
    --chrom_size_path hg19.main.nochrM.chrom.sizes   
    --obs_dim cell   
    --cpu 20   
    --chunk_size 50   
    --regions chrom100k 100000   
    --regions chrom5k 5000 
    --regions geneslop2k gencode.v28lift37.annotation.gene.bed.gz 
    --regions promoter promoter_slop2k.sorted.bed.gz 
    --regions CGI hg19_CGI.sorted.bed.gz 
    --quantifiers chrom100k count GCHN,WCGN 
    --quantifiers chrom5k count GCHN,WCGN 
    --quantifiers geneslop2k count GCHN,WCGN 
    --quantifiers promoter count GCHN,WCGN 
    --quantifiers CGI count GCHN,WCGN 
    --quantifiers chrom5k hypo-score WCGN cutoff=0.9 
    --quantifiers promoter hypo-score WCGN cutoff=0.9 
    --quantifiers CGI hypo-score WCGN cutoff=0.9 
    --quantifiers chrom5k hyper-score GCHN cutoff=0.9 
    --quantifiers promoter hyper-score GCHN cutoff=0.9 
    --quantifiers CGI hyper-score GCHN cutoff=0.9

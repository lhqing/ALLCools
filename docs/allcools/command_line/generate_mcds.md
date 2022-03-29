# allcools mcds

## Parameters 

|Required arguments                          |Descriptions                         |
|-------------------------------|-----------------------------|
|`--allc_table ALLC_TABLE`           |Contain all the ALLC file information in two tab-separated columns: 1. file_uid, 2. file_path. No header (default: None).           |
|`--output_prefix OUTPUT_PREFIX`            |Output prefix of the MCDS (default: None).            |
|`--chrom_size_path CHROM_SIZE_PATH`|Path to UCSC chrom size file. This can be generated from the genome fasta or downloaded via UCSC fetchChromSizes tools. All ALLCools functions will refer to this file whenever possible to check for chromosome names and lengths, so it is crucial to use a chrom size file consistent to the reference fasta file ever since mapping. ALLCools functions will not change or infer chromosome names. (default: None). | 
|`--mc_contexts MC_CONTEXTS [MC_CONTEXTS ...]` | Space separated mC context patterns to extract from ALLC. The context length should be the same as ALLC file context. Context pattern follows IUPAC nucleotide code, e.g. N for ATCG, H for ATC, Y for CT. (default: None).  |

<br/>

|Optional arguments                          |Descriptions                         |
|-------------------------------|-----------------------------|
|`-h, --help`           | show this help message and exit.           |
|`--split_strand`            |If true, Watson (+) and Crick (-) strands will be count separately (default: False).           |
|`--bin_sizes BIN_SIZES [BIN_SIZES ...]`|Arbitrary genomic regions can be defined in several BED files to count on. Space separated paths to each BED files, The fourth column of the BED file should be unique id of the regions (default: None). | 
|`--region_bed_paths REGION_BED_PATHS [REGION_BED_PATHS ...]` | Space separated mC context patterns to extract from ALLC. The context length should be the same as ALLC file context. Context pattern follows IUPAC nucleotide code, e.g. N for ATCG, H for ATC, Y for CT (default: None).  |
|`--region_bed_names REGION_BED_NAMES [REGION_BED_NAMES ...]` | Space separated names for each BED file provided in region_bed_paths (default: None). |
|`--cov_cutoff COV_CUTOFF` | Max cov filter for a single site in ALLC. Sites with cov > cov_cutoff will be skipped (default: 9999).  |
|`--cpu CPU` | Number of processes to use in parallel (default: 1). |
|`--max_per_mcds MAX_PER_MCDS` | Maximum number of ALLC files to aggregate into 1 MCDS, if number of ALLC provided > max_per_mcds, will generate MCDS in chunks, with same prefix provided (default: 3072). |
|`--cell_chunk_size CELL_CHUNK_SIZE` | Size of cell chunk in parallel aggregation. Do not have any effect on results. Large chunksize needs large memory (default: 100). |
|`--dtype {uint8,uint16,uint32,uint64,int8,int16,int32,int64,bool}` | Data type of MCDS count matrix. Default is np.uint32. For single cell feature count, this can be set to np.uint16 [0, 65536] to decrease file size. The values exceed min/max will be clipped while keep the mc/cov same, and a warning will be sent (default: uint32). |

## Tutorial 

**Example:**

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
GSM1922084_A02  /gale/netapp/cemba3c/projects/SingleCellMethylomes/Angermueller2016NM/allc/GSM1922084_A02.allc.tsv.gz
GSM1922085_A03  /gale/netapp/cemba3c/projects/SingleCellMethylomes/Angermueller2016NM/allc/GSM1922085_A03.allc.tsv.gz

The first part indicates the cell name (e.g. GSM1922083_A01) whereas the second part indicates the allc file path of the cell (e.g. /gale/netapp/cemba3c/projects/SingleCellMethylomes/Angermueller2016NM/allc/GSM1922083_A01.allc.tsv.gz). Make sure the two parts are separated by a space. 

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




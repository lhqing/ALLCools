from textwrap import dedent

idx_doc = (
    "If true, save an methylpy chromosome index for back compatibility. "
    "If you only use methylpy to call DMR, this don't need to be True."
)

allc_path_doc = "Path to 1 ALLC file"

allc_paths_doc = (
    "Single ALLC path contain wildcard OR multiple space separated ALLC paths "
    "OR a file contains 1 ALLC path in each row."
)

allc_table_doc = (
    "Contain all the ALLC file information in two tab-separated columns: "
    "1. file_uid, 2. file_path. No header"
)

binarize_doc = (
    "If set, binarize each single site in each individual ALLC file. "
    "This means each cytosine will only contribute at most 1 cov and 0/1 mc, "
    "this is suitable to account for single cell ALLC R1 R2 overlap issue, "
    "Only use this on single cell ALLC, not bulk ALLC."
)

bin_sizes_doc = (
    "Fix-size genomic bins can be defined by bin_sizes and chrom_size_path. "
    "Space separated sizes of genome bins, each size will be count separately."
)

bw_bin_sizes_doc = "Bin size of the BigWig files."

chrom_size_path_doc = (
    "Path to UCSC chrom size file. "
    "This can be generated from the genome fasta or downloaded via UCSC fetchChromSizes tools. "
    "All ALLCools functions will refer to this file whenever possible to check for "
    "chromosome names and lengths, so it is crucial to use a chrom size file consistent "
    "to the reference fasta file ever since mapping. "
    "ALLCools functions will not change or infer chromosome names."
)

compress_level_doc = "Compression level for the output file"

cov_cutoff_doc = "Max cov filter for a single site in ALLC. Sites with cov > cov_cutoff will be skipped."

cpu_basic_doc = "Number of processes to use in parallel."

mc_contexts_doc = (
    "Space separated mC context patterns to extract from ALLC. "
    "The context length should be the same as ALLC file context. "
    "Context pattern follows IUPAC nucleotide code, e.g. N for ATCG, H for ATC, Y for CT."
)

mc_context_mcad_doc = (
    "mC context pattern to extract from ALLC. "
    "Context pattern follows IUPAC nucleotide code, e.g. N for ATCG, H for ATC, Y for CT."
    "Note that generate_mcad only take one mC context"
)

reference_fasta_doc = (
    "Path to 1 genome reference FASTA file (the one used for mapping), "
    "use samtools fadix to build .fai index first. Do not compress that file."
)

region_bed_names_doc = (
    "Space separated names for each BED file provided in region_bed_paths."
)

region_bed_paths_doc = (
    "Arbitrary genomic regions can be defined in several BED files to count on. "
    "Space separated paths to each BED files, "
    "The fourth column of the BED file should be unique id of the regions."
)

region_bed_path_mcad_doc = (
    "Arbitrary genomic regions can be defined in one BED file to count on. "
    "The fourth column of the BED file should be unique id of the regions."
)

region_doc = (
    "Only extract records from certain genome region(s) via tabix, "
    "multiple region can be provided in tabix form. If region is not None, will not run in parallel"
)

remove_additional_chrom_doc = (
    "Whether to remove rows with unknown chromosome instead of raising KeyError"
)

rna_table_doc = (
    "This is only for mCT data when we have RNA BAM file for each single cell. "
    "Contain all the RNA BAM file information in 2 columns: 1. file_uid, 2. file_path. No header."
)

snp_doc = "If true, means the input allc contain snp information, and the allc processing will take care that."

split_strand_doc = "If true, Watson (+) and Crick (-) strands will be count separately"

strandness_doc = (
    "What to do with strand information, possible values are: "
    "1. both: save +/- strand together in one file; "
    "2. split: save +/- strand into two separate files, with suffix contain Watson (+) and Crick (-); "
    "3. merge: This will only merge the count on adjacent CpG in +/- strands, "
    "only work for CpG like context. For non-CG context, its the same as both."
)

generate_dataset_doc = (
    "Generate MCDS dataset from a list of ALLC files (recorded in the allc_table). "
    "Multiple region sets, methylation contexts and quantification types can be included in one command."
)

generate_dataset_obs_dim_doc = 'Name of the observation dimension.'

generate_dataset_chunk_size_doc = 'Chunk allc_table with chunk_size when generate dataset in parallel'

generate_dataset_regions_doc = (
    'Definition of genomic regions in the form of "--regions {region_name} {region_definition}". '
    'This parameter can be specified multiple times, to allow quantification of multiple region sets '
    'in the same MCDS dataset. Several cases are allowed: '
    '1) a integer number means fix-sized genomic bins, region bed and region id will be generated '
    'automatically based on the chrom_size_path parameter (e.g., "--regions chrom100k 100000"); '
    '2) a path to a three-column bed file, in this case, '
    'a forth column containing region id in the form of {region_name}_{i} will be added automatically '
    '(e.g., "--regions gene /path/to/gene_bed_no_id.bed", '
    'where the bed file only has chrom, start, end columns); '
    '3) a path to a four-column bed file, in this case, the forth column will be treated as region id '
    'and the region ids must be UNIQUE. (e.g., "--regions gene /path/to/gene_bed_with_id.bed", '
    'where the bed file has chrom, start, end, id columns).'
)

generate_dataset_quantifiers_doc = (
    'Definition of genome region quantifiers in the form of '
    '"--quantifiers {region_name} {quant_type} {mc_contexts} {optional_parameter}". '
    'The region_name determines which region set this quantifier applies to, '
    'region_name must be defined by "--regions" parameter. '
    'The quant_type specify which quantifiers, it must be in ["count", "hypo-score", "hyper-score"]. '
    'The mc_contexts specify a comma separated mC context list, '
    'it must be the same size as the ALLC table, and uses IUPAC base abbreviation. '
    '--quantifiers parameter can be specified multiple times, '
    'to allow different quantification for different region sets, '
    'or multiple quantification for the same region set. '
    'Some examples: '
    '1) To quantify raw counts of a region set in mCG and mCH context: '
    '"--quantifiers gene count CGN,CHN" '
    '2) To quantify the mCG hypo-methylation score of chrom 5Kb bins: '
    '"--quantifiers chrom5k hypo-score CGN cutoff=0.9", '
    'by default, cutoff=0.9, so the last part is optional. '
    '3) To ALSO quantify the mCG raw counts of chrom 5Kb bins in the same MCDS, '
    'just specify another quantifiers in the same command: '
    '"--quantifiers chrom5k count CGN", note the count matrix of chrom5k will be large. '
    'Its not usually needed, but you have the option if needed.'
)

table_to_allc_doc = 'Convert different kinds of methylation table into ALLC format. ' \
                    'Currently, only plain text table is accepted.'
table_to_allc_input_path = 'input path of the table'
table_to_allc_output_prefix = 'output prefix of the ALLC table'
table_to_allc_sep = 'character to separate columns in the table'
table_to_allc_header = 'Whether the table contains header line or not'
table_to_allc_chunk_size = 'chunk_size to perform conversion'
table_to_allc_chrom = 'the chromosome column number, 0-based index'
table_to_allc_pos = 'the position column number, 0-based index'
table_to_allc_strand = 'the strand column number, 0-based index. ' \
                       'If not provided, will infer automatically based on the fasta_path'
table_to_allc_context = 'the cytosine context column number, 0-based index. ' \
                        'If not provided, will inter automatically based on the fasta_path'
table_to_allc_mc = 'the methylated cytosine count column number, 0-based index.'
table_to_allc_uc = 'the unmethylated cytosine count column number, 0-based index.'
table_to_allc_cov = 'the total cytosine coverage count column number, 0-based index.'
table_to_allc_mc_frac = 'the methylation fraction column number, 0-based index.'
table_to_allc_pseudo_count = 'Use this pseudo_count number as the total cytosine coverage count, ' \
                             'if the "cov" column is missing and "mc_frac" column is provided.'
table_to_allc_fasta_path = 'the genome FASTA file path, ' \
                           'required if either "strand" or "context" column is missing.'
table_to_allc_num_upstream_bases = 'number of up stream bases to include when get cytosine context.'
table_to_allc_num_downstream_bases = 'number of down stream bases to include when get cytosine context.'
table_to_allc_add_chr = 'whether add "chr" before the chromosome name.'
table_to_allc_sort = 'whether sort the ALLC table after conversion.'


def doc_params(**kwds):
    """\
    Docstrings should start with "\" in the first line for proper formatting.
    """

    def dec(obj):
        obj.__doc__ = dedent(obj.__doc__).format(**kwds)
        return obj

    return dec

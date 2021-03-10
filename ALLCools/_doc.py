from textwrap import dedent

idx_doc = "If true, save an methylpy chromosome index for back compatibility. " \
          "If you only use methylpy to call DMR, this don't need to be True."

allc_path_doc = "Path to 1 ALLC file"

allc_paths_doc = "Single ALLC path contain wildcard OR multiple space separated ALLC paths " \
                 "OR a file contains 1 ALLC path in each row."

allc_table_doc = "Contain all the ALLC file information in two tab-separated columns: " \
                 "1. file_uid, 2. file_path. No header"

binarize_doc = 'If set, binarize each single site in each individual ALLC file. ' \
               'This means each cytosine will only contribute at most 1 cov and 0/1 mc, ' \
               'this is suitable to account for single cell ALLC R1 R2 overlap issue, ' \
               'Only use this on single cell ALLC, not bulk ALLC.'

bin_sizes_doc = "Fix-size genomic bins can be defined by bin_sizes and chrom_size_path. " \
                "Space separated sizes of genome bins, each size will be count separately."

chrom_size_path_doc = "Path to UCSC chrom size file. " \
                      "This can be generated from the genome fasta or downloaded via UCSC fetchChromSizes tools. " \
                      "All ALLCools functions will refer to this file whenever possible to check for " \
                      "chromosome names and lengths, so it is crucial to use a chrom size file consistent " \
                      "to the reference fasta file ever since mapping. " \
                      "ALLCools functions will not change or infer chromosome names."

compress_level_doc = "Compression level for the output file"

cov_cutoff_doc = "Max cov filter for a single site in ALLC. Sites with cov > cov_cutoff will be skipped."

cpu_basic_doc = 'Number of processes to use in parallel.'

mc_contexts_doc = "Space separated mC context patterns to extract from ALLC. " \
                  "The context length should be the same as ALLC file context. " \
                  "Context pattern follows IUPAC nucleotide code, e.g. N for ATCG, H for ATC, Y for CT."

mc_context_mcad_doc = 'mC context pattern to extract from ALLC. ' \
                      'Context pattern follows IUPAC nucleotide code, e.g. N for ATCG, H for ATC, Y for CT.'

reference_fasta_doc = 'Path to 1 genome reference FASTA file (the one used for mapping), ' \
                      'use samtools fadix to build .fai index first. Do not compress that file.'

region_bed_names_doc = 'Space separated names for each BED file provided in region_bed_paths.'

region_bed_paths_doc = 'Arbitrary genomic regions can be defined in several BED files to count on. ' \
                       'Space separated paths to each BED files, ' \
                       'The fourth column of the BED file should be unique id of the regions.'

region_bed_path_mcad_doc = 'Arbitrary genomic regions can be defined in one BED file to count on. ' \
                           'The fourth column of the BED file should be unique id of the regions.'

region_doc = "Only extract records from certain genome region(s) via tabix, " \
             "multiple region can be provided in tabix form. If region is not None, will not run in parallel"

remove_additional_chrom_doc = "Whether to remove rows with unknown chromosome instead of raising KeyError"

rna_table_doc = "This is only for mCT data when we have RNA BAM file for each single cell. " \
                "Contain all the RNA BAM file information in 2 columns: 1. file_uid, 2. file_path. No header."

snp_doc = 'If true, means the input allc contain snp information, and the allc processing will take care that.'

split_strand_doc = 'If true, Watson (+) and Crick (-) strands will be count separately'

strandness_doc = "What to do with strand information, possible values are: " \
                 "1. both: save +/- strand together in one file; " \
                 "2. split: save +/- strand into two separate files, with suffix contain Watson (+) and Crick (-); " \
                 "3. merge: This will only merge the count on adjacent CpG in +/- strands, " \
                 "only work for CpG like context. For non-CG context, its the same as both."


def doc_params(**kwds):
    """\
    Docstrings should start with "\" in the first line for proper formatting.
    """

    def dec(obj):
        obj.__doc__ = dedent(obj.__doc__).format(**kwds)
        return obj

    return dec

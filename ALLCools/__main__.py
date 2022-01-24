"""
CLI defined here

When adding new function:
1. add a func_register_subparser function to register the subparser
2. add a condition in main func about this new func name, import the real func as func in main
"""

import argparse
import inspect
import logging
import subprocess
import sys

import ALLCools
from ._doc import *
import os

os.environ["NUMEXPR_MAX_THREADS"] = "8"

log = logging.getLogger()

DESCRIPTION = """
The ALLCools command line toolkit contains multiple functions to manipulate the ALLC format, 
a core file format that stores single base level methylation information.
Throughout this toolkit, we use bgzip/tabix to compress and index the ALLC file to allow 
flexible data query from the ALLC file.

Current Tool List in ALLCools:

[Generate ALLC]
bam-to-allc          - Generate 1 ALLC file from 1 position sorted BAM file via 
                       samtools mpileup.

[Manipulate ALLC]
standardize-allc     - Validate 1 ALLC file format, standardize the chromosome names, 
                       compression format (bgzip) and index (tabix).
tabix-allc           - A simple wrapper of tabix command to index 1 ALLC file.
profile-allc         - Generate some summary statistics of 1 ALLC
merge-allc           - Merge N ALLC files into 1 ALLC file
extract-allc         - Extract information (strand, context) from 1 ALLC file

[Get Region Level]
allc-to-bigwig       - Generate coverage (cov) and ratio (mc/cov) bigwig track files 
                       from 1 ALLC file
allc-to-region-count - Count region level mc, cov by genome bins or provided BED files.
generate-mcds        - Generate methylation dataset (MCDS) for a group of ALLC file and 
                       different region sets. This is a convenient wrapper function for 
                       a bunch of allc-to-region-count and xarray integration codes. 
                       MCDS is inherit from xarray.DataSet
generate-mcad        - Generate mCG hypo-methylation score AnnData dataset (MCAD) for 
                       a group of ALLC file and one region set.
"""

EPILOG = """
Author: Hanqing Liu

See ALLCools documentation here: https://lhqing.github.io/ALLCools/intro.html
"""


class NiceFormatter(logging.Formatter):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).
    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = "{}: {}".format(record.levelname, record.msg)
        return super().format(record)


def validate_environment():
    try:
        subprocess.run(
            ["tabix", "--version"],
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            encoding="utf8",
            check=True,
        )
        subprocess.run(
            ["bgzip", "--version"],
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            encoding="utf8",
            check=True,
        )
        subprocess.run(
            ["bedtools", "--version"],
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            encoding="utf8",
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(e.stderr)
        raise e
    return True


def setup_logging(stdout=False, quiet=False, debug=False):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Attach handler to the global logger object
    """
    # Due to backwards compatibility, logging output is sent to standard output
    # instead of standard error if the -o option is used.
    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setFormatter(NiceFormatter())
    # debug overrides quiet
    if debug:
        level = logging.DEBUG
    elif quiet:
        level = logging.ERROR
    else:
        level = logging.INFO
    stream_handler.setLevel(level)
    log.setLevel(level)
    log.addHandler(stream_handler)


def _str_to_bool(v: str) -> bool:
    if v.lower() in {
        "true",
        "1",
        "t",
        "y",
        "yes",
        "yeah",
        "yup",
        "certainly",
        "uh-huh",
    }:
        return True
    else:
        return False


def bam_to_allc_register_subparser(subparser):
    parser = subparser.add_parser(
        "bam-to-allc",
        aliases=["allc", "2allc"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Take 1 position sorted BAM file, generate 1 ALLC file.",
    )

    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument(
        "--bam_path",
        "-bam",
        dest="bam_path",
        type=str,
        required=True,
        help="Path to 1 position sorted BAM file",
    )

    parser_req.add_argument(
        "--reference_fasta", type=str, required=True, help=reference_fasta_doc
    )

    parser_req.add_argument(
        "--output_path", type=str, required=True, help="Path to output ALLC file"
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help="Number of processes to use in parallel. DO NOT use cpu > 1 for single cell ALLC generation. "
             "Parallel on cell level is better for single cell project.",
    )

    parser.add_argument(
        "--num_upstr_bases",
        type=int,
        default=0,
        help="Number of upstream base(s) of the C base to include in ALLC context column, "
             "usually use 0 for normal BS-seq, 1 for NOMe-seq.",
    )

    parser.add_argument(
        "--num_downstr_bases",
        type=int,
        default=2,
        help="Number of downstream base(s) of the C base to include in ALLC context column, "
             "usually use 2 for both BS-seq and NOMe-seq.",
    )

    parser.add_argument(
        "--min_mapq",
        type=int,
        default=10,
        help="Minimum MAPQ for a read being considered, samtools mpileup parameter, see samtools documentation.",
    )

    parser.add_argument(
        "--min_base_quality",
        type=int,
        default=20,
        help="Minimum base quality for a base being considered, samtools mpileup parameter, "
             "see samtools documentation.",
    )

    parser.add_argument(
        "--compress_level", type=int, default=5, help=compress_level_doc
    )

    parser.add_argument(
        "--save_count_df",
        dest="save_count_df",
        action="store_true",
        help="If present, save an ALLC context count table next to ALLC file.",
    )
    parser.set_defaults(save_count_df=False)


def standardize_allc_register_subparser(subparser):
    parser = subparser.add_parser(
        "standardize-allc",
        aliases=["standard"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Standardize 1 ALLC file by checking: "
             "1. No header in the ALLC file; "
             "2. Chromosome names in ALLC must be exactly same as those "
             "in the chrom_size_path file; "
             "3. Output file will be bgzipped with .tbi index; "
             "4. Remove additional chromosome (remove_additional_chrom=True) or "
             "raise KeyError if unknown chromosome found (default)",
    )
    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument("--allc_path", type=str, required=True, help=allc_path_doc)

    parser.add_argument(
        "--chrom_size_path", type=str, required=True, help=chrom_size_path_doc
    )

    parser.add_argument(
        "--compress_level", type=int, required=False, default=5, help=compress_level_doc
    )

    parser.add_argument(
        "--remove_additional_chrom",
        dest="remove_additional_chrom",
        action="store_true",
        help=remove_additional_chrom_doc,
    )
    parser.set_defaults(remove_additional_chrom=False)


def tabix_allc_register_subparser(subparser):
    parser = subparser.add_parser(
        "tabix-allc",
        aliases=["tbi"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="a simple wrapper of tabix command to index 1 ALLC file",
    )
    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument(dest="allc_path", nargs="?", type=str, help=allc_path_doc)

    parser.add_argument(
        "--reindex",
        dest="reindex",
        action="store_true",
        help="If present, will force regenerate the ALLC index.",
    )
    parser.set_defaults(reindex=False)


def profile_allc_register_subparser(subparser):
    parser = subparser.add_parser(
        "profile-allc",
        aliases=["profile"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Generate some summary statistics of 1 ALLC.",
    )

    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument("--allc_path", type=str, required=True, help=allc_path_doc)

    parser.add_argument(
        "--keep_n",
        dest="drop_n",
        action="store_false",
        help="Whether to keep context that contain N, such as CCN. "
             "This is usually very rare and need to be dropped.",
    )
    parser.set_defaults(drop_n=True)

    parser.add_argument(
        "--n_rows",
        type=int,
        default=1000000,
        help="Number of rows to calculate the profile from. "
             "The default number is usually sufficient to get pretty precise assumption.",
    )

    parser.add_argument(
        "--output_path",
        type=str,
        default=None,
        help="Path of the output file. If None, will save the profile next to input ALLC file.",
    )


def merge_allc_register_subparser(subparser):
    parser = subparser.add_parser(
        "merge-allc",
        aliases=["merge"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Merge N ALLC files into 1 ALLC file",
    )

    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument(
        "--allc_paths", type=str, nargs="+", required=True, help=allc_paths_doc
    )

    parser_req.add_argument(
        "--output_path",
        type=str,
        required=True,
        help="Path to the output merged ALLC file.",
    )

    parser_req.add_argument(
        "--chrom_size_path", type=str, required=True, help=chrom_size_path_doc
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=10,
        help=f"{cpu_basic_doc} The real CPU usage is ~1.5 times than this number, "
             f"due to the sub processes of handling ALLC files using tabix/bgzip. "
             f"Monitor the CPU and Memory usage when running this function.",
    )

    parser.add_argument(
        "--bin_length",
        type=int,
        default=10000000,
        help="Length of the genome bin in each parallel job, large number means more memory usage.",
    )

    parser.add_argument("--snp", dest="snp", action="store_true", help=snp_doc)
    parser.set_defaults(snp=False)


def extract_context_allc_register_subparser(subparser):
    parser = subparser.add_parser(
        "extract-allc",
        aliases=["extract"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Extract information (strand, context) from 1 ALLC file. "
             "Able to save to several different format.",
    )

    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument("--allc_path", type=str, required=True, help=allc_path_doc)

    parser_req.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="Path prefix of the output ALLC file.",
    )

    parser_req.add_argument(
        "--mc_contexts", type=str, required=True, nargs="+", help=mc_contexts_doc
    )

    parser_req.add_argument(
        "--chrom_size_path", type=str, required=True, help=chrom_size_path_doc
    )

    parser.add_argument(
        "--strandness",
        type=str,
        default="both",
        choices=["both", "split", "merge"],
        help="What to do with strand information, possible values are: "
             "1. both: save +/- strand together in one file without any modification; "
             "2. split: save +/- strand into two separate files, with suffix contain Watson (+) and Crick (-); "
             "3. merge: This will only merge the count on adjacent CpG in +/- strands, only work for CpG like context. "
             "For non-CG context, its the same as both.",
    )

    parser.add_argument(
        "--output_format",
        type=str,
        default="allc",
        choices=["allc", "bed5"],
        help="Output format of extracted information, possible values are: "
             "1. allc: keep the allc format; "
             "2. bed5: 5-column bed format, chrom, pos, pos, mc, cov; ",
    )

    parser.add_argument("--region", type=str, default=None, help=region_doc)

    parser.add_argument("--cov_cutoff", type=int, default=99999, help=cov_cutoff_doc)

    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help=cpu_basic_doc + " This function parallel on region level and "
                             "will generate a bunch of small files if cpu > 1. "
                             "Do not use cpu > 1 for single cell region count. "
                             "For single cell data, parallel on cell level is better.",
    )
    return


def allc_to_region_count_register_subparser(subparser):
    parser = subparser.add_parser(
        "allc-to-region-count",
        aliases=["region", "2region"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Calculate mC and cov at regional level. Region can be provided in 2 forms: "
             "1. BED file, provided by region_bed_paths, "
             "containing arbitrary regions and use bedtools map to calculate; "
             "2. Fix-size non-overlap genome bins, provided by bin_sizes, "
             "Form 2 is much faster to calculate than form 1. "
             "The output file is in 6-column bed-like format: "
             "chrom start end region_uid mc cov",
    )

    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument(
        "--allc_path",
        "-allc",
        dest="allc_path",
        type=str,
        required=True,
        help=allc_path_doc,
    )

    parser_req.add_argument(
        "--output_prefix",
        "-out",
        dest="output_prefix",
        type=str,
        required=True,
        help="Path prefix of the output region count file.",
    )

    parser_req.add_argument(
        "--chrom_size_path", type=str, required=True, help=chrom_size_path_doc
    )

    parser_req.add_argument(
        "--mc_contexts", type=str, nargs="+", required=True, help=mc_contexts_doc
    )

    parser.add_argument(
        "--split_strand",
        dest="split_strand",
        action="store_true",
        help=split_strand_doc,
    )
    parser.set_defaults(split_strand=False)

    parser.add_argument(
        "--region_bed_paths",
        type=str,
        nargs="+",
        default=None,
        help=region_bed_paths_doc,
    )

    parser.add_argument(
        "--region_bed_names",
        type=str,
        nargs="+",
        default=None,
        help=region_bed_names_doc,
    )

    parser.add_argument(
        "--bin_sizes", type=int, nargs="+", default=None, help=bin_sizes_doc
    )

    parser.add_argument("--cov_cutoff", type=int, default=9999, help=cov_cutoff_doc)

    parser.add_argument(
        "--save_zero_cov",
        dest="save_zero_cov",
        action="store_true",
        help="If present, save the regions that have 0 cov in output, "
             "only apply to region count but not the chromosome count.",
    )
    parser.set_defaults(save_zero_cov=False)

    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help=cpu_basic_doc + " This function parallel on region level and "
                             "will generate a bunch of small files if cpu > 1. "
                             "Do not use cpu > 1 for single cell region count. "
                             "For single cell data, parallel on cell level is better.",
    )
    return


def allc_to_bigwig_register_subparser(subparser):
    parser = subparser.add_parser(
        "allc-to-bigwig",
        aliases=["bw", "2bw"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Generate bigwig file(s) from 1 ALLC file.",
    )

    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument("--allc_path", type=str, required=True, help=allc_path_doc)

    parser_req.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="Path prefix of the output BigWig file.",
    )

    parser.add_argument("--bin_size", type=int, default=50, help=bw_bin_sizes_doc)

    parser_req.add_argument(
        "--mc_contexts", type=str, nargs="+", required=True, help=mc_contexts_doc
    )

    parser_req.add_argument(
        "--chrom_size_path", type=str, required=True, help=chrom_size_path_doc
    )

    parser.add_argument(
        "--strandness",
        type=str,
        choices=["split", "both"],
        default="both",
        help=strandness_doc,
    )
    return


def generate_mcds_register_subparser(subparser):
    parser = subparser.add_parser(
        "generate-mcds",
        aliases=["mcds"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Generate MCDS from ALLC files and region sets.",
    )

    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument(
        "--allc_table", type=str, required=True, help=allc_table_doc
    )

    parser_req.add_argument(
        "--output_prefix", type=str, required=True, help="Output prefix of the MCDS"
    )

    parser_req.add_argument(
        "--chrom_size_path", type=str, required=True, help=chrom_size_path_doc
    )

    parser_req.add_argument(
        "--mc_contexts", type=str, required=True, nargs="+", help=mc_contexts_doc
    )

    parser.add_argument(
        "--split_strand",
        dest="split_strand",
        action="store_true",
        help=split_strand_doc,
    )
    parser.set_defaults(split_strand=False)

    parser.add_argument(
        "--bin_sizes", type=int, nargs="+", default=None, help=bin_sizes_doc
    )

    parser.add_argument(
        "--region_bed_paths",
        type=str,
        nargs="+",
        default=None,
        help=region_bed_paths_doc,
    )

    parser.add_argument(
        "--region_bed_names",
        type=str,
        nargs="+",
        default=None,
        help=region_bed_names_doc,
    )

    parser.add_argument("--cov_cutoff", type=int, default=9999, help=cov_cutoff_doc)

    parser.add_argument("--cpu", type=int, default=1, help=cpu_basic_doc)

    parser.add_argument(
        "--max_per_mcds",
        type=int,
        default=3072,
        help="Maximum number of ALLC files to aggregate into 1 MCDS, if number of ALLC provided > max_per_mcds, "
             "will generate MCDS in chunks, with same prefix provided.",
    )

    parser.add_argument(
        "--cell_chunk_size",
        type=int,
        default=100,
        help="Size of cell chunk in parallel aggregation. Do not have any effect on results. "
             "Large chunksize needs large memory.",
    )

    parser.add_argument(
        "--dtype",
        type=str,
        default="uint32",
        choices=[
            "uint8",
            "uint16",
            "uint32",
            "uint64",
            "int8",
            "int16",
            "int32",
            "int64",
            "bool",
        ],
        help="Data type of MCDS count matrix. Default is np.uint32. "
             "For single cell feature count, this can be set to np.uint16 [0, 65536] to decrease file size. "
             "The values exceed min/max will be clipped while keep the mc/cov same, and a warning will be sent.",
    )

    parser.add_argument(
        "--engine",
        type=str,
        default="zarr",
        choices=["zarr", "netcdf"],
        help="use zarr or netcdf to store dataset, default is zarr"
    )
    return


def generate_mcad_register_subparser(subparser):
    parser = subparser.add_parser(
        "generate-mcad",
        aliases=["mcad"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Generate MCAD from ALLC files and one region set.",
    )

    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument(
        "--allc_table", type=str, required=True, help=allc_table_doc
    )

    parser_req.add_argument(
        "--bed_path", type=str, required=True, help=region_bed_path_mcad_doc
    )

    parser_req.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help='Output prefix of the MCAD, a suffix ".mcad" will be added.',
    )

    parser_req.add_argument(
        "--mc_context", type=str, required=True, help=mc_context_mcad_doc
    )

    parser.add_argument("--cpu", type=int, default=1, help=cpu_basic_doc)

    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.9,
        help="Values smaller than cutoff will be stored as 0, which reduces the file size.",
    )

    parser.add_argument(
        "--reverse_value",
        dest="reverse_value",
        action="store_true",
        help="If present, use cdf instead of sf to make hyper-methylation events having higher values.",
    )
    parser.set_defaults(reverse_value=False)
    return


def convert_mcds_to_zarr_register_subparser(subparser):
    parser = subparser.add_parser(
        "convert-mcds-to-zarr",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Convert xarray.Dataset stored in other backends into zarr backend.",
    )

    parser.add_argument(
        "paths",
        type=str,
        nargs='+',
        help='provide one or multiple MCDS paths, wildcards supported'
    )
    return


def generate_dataset_register_subparser(subparser):
    parser = subparser.add_parser(
        "generate-dataset",
        aliases=["dataset"],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="",
    )

    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument(
        "--allc_table",
        type=str,
        required=True,
        help=''
    )

    parser_req.add_argument(
        "--output_path",
        type=str,
        required=True,
        help=''
    )

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True,
        help=''
    )

    parser.add_argument(
        "--obs_dim",
        type=str,
        default='cell',
        help=''
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help=''
    )

    parser.add_argument(
        "--chunk_size",
        type=int,
        default=10,
        help=''
    )

    parser.add_argument(
        "--regions",
        type=str,
        nargs=2,
        action='append',
        help=''
    )

    parser.add_argument(
        "--quantifiers",
        type=str,
        nargs='+',
        action='append',
        help=''
    )
    return


def main():
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(title="functions", dest="command", metavar="")

    # add subparsers
    current_module = sys.modules[__name__]
    # get all functions in parser
    for name, register_subparser_func in inspect.getmembers(
            current_module, inspect.isfunction
    ):
        if "register_subparser" in name:
            register_subparser_func(subparsers)

    # initiate
    args = None
    if len(sys.argv) > 1:
        # print out version
        if sys.argv[1] in ["-v", "--version"]:
            print(ALLCools.__version__)
            exit()
        else:
            args = parser.parse_args()
    else:
        # print out help
        parser.parse_args(["-h"])
        exit()

    # set up logging
    if not logging.root.handlers:
        setup_logging(stdout=True, quiet=False)

    # execute command
    args_vars = vars(args)
    for k, v in args_vars.items():
        log.info(f"{k}\t{v}")

    cur_command = args_vars.pop("command").lower().replace("_", "-")
    # Do real import here:
    if cur_command in ["bam-to-allc", "allc", "2allc"]:
        from ._bam_to_allc import bam_to_allc as func
    elif cur_command in ["standardize-allc", "standard"]:
        from .utilities import standardize_allc as func
    elif cur_command in ["tabix-allc", "tbi"]:
        from .utilities import tabix_allc as func
    elif cur_command in ["profile-allc", "profile"]:
        from .utilities import profile_allc as func
    elif cur_command in ["merge-allc", "merge"]:
        from ._merge_allc import merge_allc_files as func
    elif cur_command in ["extract-allc", "extract"]:
        from ._extract_allc import extract_allc as func
    elif cur_command in ["allc-to-region-count", "region", "2region"]:
        from ._allc_to_region_count import allc_to_region_count as func
    elif cur_command in ["allc-to-bigwig", "bw", "2bw"]:
        from ._allc_to_bigwig_new import allc_to_bigwig as func
        # from ._allc_to_bigwig import allc_to_bigwig as func  # old version using ucsc tools
    elif cur_command in ["generate-mcds", "mcds"]:
        from .count_matrix.mcds import generate_mcds as func
    elif cur_command in ['convert-mcds-to-zarr']:
        from .mcds.utilities import convert_to_zarr as func
    elif cur_command in ["generate-mcad", "mcad"]:
        from .count_matrix.mcad import generate_mcad as func
    elif cur_command in ["generate-dataset", "dataset"]:
        from .count_matrix.dataset import generate_dataset as func
    else:
        log.debug(f"{cur_command} is not an valid sub-command")
        parser.parse_args(["-h"])
        return

    validate_environment()

    # run the command
    log.info(f"Executing {cur_command}...")
    func(**args_vars)
    log.info(f"{cur_command} finished.")
    return


if __name__ == "__main__":
    main()

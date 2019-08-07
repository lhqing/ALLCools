import pathlib
from .utilities import get_fasta

def parallel_ame(bed_files, background_files, motif_files,
                 output_dir, reference_fasta, cpu=1, slop_b=50):
    """
    Run motif enrichment analysis on EACH bed_file with ALL the motif_files,
    using ALL background_files (merged) as control region

    Parameters
    ----------
    bed_files
    background_files
    motif_files
    cpu

    Returns
    -------

    """
    output_dir = pathlib.Path(output_dir)
    if isinstance(bed_files, str):
        bed_files = [bed_files]
    if isinstance(background_files, str):
        background_files = [background_files]

    # get bed FASTQ
    fasta_list = []
    for bed_file in bed_files:
        bed_file = pathlib.Path(bed_file)
        output_fasta_path = output_dir / (bed_file.name + '.fasta')
        get_fasta([bed_file],
                  fasta_path=reference_fasta,
                  output_path=output_fasta_path,
                  slop_b=slop_b,
                  chrom_size_path=None,
                  cpu=cpu, sort_mem_gbs=min(int(cpu * 3), 30))
        fasta_list.append(output_fasta_path)

    # get merged background FASTQ
    control_fasta = output_dir / 'CONTROL_REGIONS.fasta'
    get_fasta(background_files,
              fasta_path=reference_fasta,
              output_path=control_fasta,
              slop_b=slop_b,
              chrom_size_path=None,
              cpu=cpu, sort_mem_gbs=min(int(cpu * 3), 30))

    # split motif files based on cpu


    # run single_bed_single_motif_runner

    # merge results

    return


def single_bed_single_motif_runner():
    return

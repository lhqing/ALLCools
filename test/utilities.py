import os
import pathlib
import shlex
import subprocess


def data_file_path(path):
    return os.path.join(os.path.dirname(__file__), 'file', path)


def output_path(path):
    return os.path.join(os.path.dirname(__file__), 'output', path)


def run_command(cmd):
    try:
        subprocess.run(shlex.split(cmd),
                       check=True,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, encoding='utf8')
    except subprocess.CalledProcessError as e:
        print(e.stdout)
        print(e.stderr)
        raise e


# mkdir if not exist
pathlib.Path(os.path.join(os.path.dirname(__file__), 'output')).mkdir(exist_ok=True)

REFERENCE_FASTA = '/home/hanliu/ref/mouse/genome/fasta/with_chrl/mm10_with_chrl.fa'
CHROM_SIZE_PATH = data_file_path('mm10.main.chrom.sizes')
GENE_BED = '/home/hanliu/ref/mouse/gencode/vm16/gencode.vM16.annotation.gene.bed.gz'

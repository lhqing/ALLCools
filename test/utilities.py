import os
import shlex
import subprocess

REFERENCE_FASTA = '/gale/netapp/home/hanliu/ref/mouse/genome/fasta/with_chrl/mm10_with_chrl.fa'


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

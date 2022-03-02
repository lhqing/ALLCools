import pathlib

import xarray as xr
from Bio import SeqIO


def prepare_motif_scan_snakemake(output_dir,
                                 fasta_path,
                                 region_dim,
                                 motif_dim='motif',
                                 motif_set_path=None,
                                 chunk_size=100000,
                                 combine_cluster=True,
                                 fnr_fpr_fold=1000,
                                 cpu=10):
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # split FASTA chunks
    with open(fasta_path) as fasta:
        seqs = []
        for i, seq in enumerate(SeqIO.parse(fasta, "fasta")):
            seqs.append(seq)
            if (i + 1) % chunk_size == 0:
                motif_chunk_dir = output_dir / f'chunk_{i // chunk_size:06d}'
                motif_chunk_dir.mkdir(exist_ok=True)
                with open(motif_chunk_dir / f'regions.fasta', 'w') as out_f:
                    SeqIO.write(seqs, out_f, 'fasta')
                seqs = []

        if len(seqs) > 0:
            motif_chunk_dir = output_dir / f'chunk_{i // chunk_size:06d}'
            motif_chunk_dir.mkdir(exist_ok=True)
            with open(motif_chunk_dir / f'regions.fasta', 'w') as out_f:
                SeqIO.write(seqs, out_f, 'fasta')

    # make snakefile
    snakemake_cmds = []
    for motif_chunk_dir in output_dir.glob('chunk_*'):
        snakemake_string = f"""
region_dim = '{region_dim}'
cpu = {cpu}
combine_cluster = {combine_cluster}
dim = '{motif_dim}'
motif_set_path = {motif_set_path}
fnr_fpr_fold = {fnr_fpr_fold}


from ALLCools.motif import MotifSet

if motif_set_path is not None:
    motif_set = joblib.load(motif_set_path)
else:
    # default motif database from three sources
    from ALLCools.motif import get_default_motif_set
    motif_set = get_default_motif_set()

if fnr_fpr_fold != 1000:
    motif_set.calculate_threshold(
        cpu=cpu, method="balance", threshold_value=fnr_fpr_fold
    )

rule run_scan:
    input:
        fasta_path='regions.fasta',
    output:
        motif_dir=directory('{region_dim}_{motif_dim}/'),
        motif_cluster_dir=directory('{region_dim}_{motif_dim}-cluster/'),
    run:
        motif_ds = motif_set.scan_motifs(
            fasta_path='regions.fasta',
            output_dir='./',
            cpu=cpu,
            region_dim='{region_dim}',
            combine_cluster=combine_cluster,
            motif_dim='{motif_dim}',
            chunk_size=120,
        )
        with open('Success', 'w') as f:
            f.write('42')
"""
        with open(motif_chunk_dir / 'Snakefile', 'w') as f:
            f.write(snakemake_string)
            motif_chunk_dir = motif_chunk_dir.absolute()
            snakemake_cmd = f'snakemake -d {motif_chunk_dir} -s {motif_chunk_dir}/Snakefile -j 1'
            snakemake_cmds.append(snakemake_cmd)

    with open(f'{output_dir}/snakemake_cmds.txt', 'w') as f:
        f.write('\n'.join(snakemake_cmds))

    with open(f'{output_dir}/status', 'w') as f:
        f.write('prepared')

    print('Snakemake files generated for motif scan.')
    print(f'Snakemake commands are listed in\n{output_dir}/snakemake_cmds.txt')
    print('Execute the snakemake commands via proper scheduler or work station.')
    return


def check_snakemake_success(output_dir):
    output_dir = pathlib.Path(output_dir)
    all_success = True
    na_paths = []
    for chunk_dir in output_dir.glob('chunk_*'):
        success_flag = chunk_dir / 'Success'
        if not success_flag.exists():
            all_success = False
            na_paths.append(str(chunk_dir))
    if all_success:
        return True
    else:
        path_str = "\n".join(na_paths)
        print(f'These chunks have not been executed successfully:')
        print(path_str)
        return False


def save_motif_chunks(motif_chunk_dir, region_dim, output_path, is_motif_cluster):
    if is_motif_cluster:
        chunk_list = sorted(pathlib.Path(motif_chunk_dir).glob('chunk_*/dmr_motif-cluster'))
    else:
        chunk_list = sorted(pathlib.Path(motif_chunk_dir).glob('chunk_*/dmr_motif'))

    # load all chunks, determine coordinate dtypes
    total_ds = xr.open_mfdataset(chunk_list,
                                 engine='zarr',
                                 concat_dim=region_dim,
                                 combine='nested')
    coord_dtypes = {
        # specifically deal with string dtype
        k: v.dtype if v.dtype != 'O' else v.astype(str).dtype
        for k, v in total_ds.coords.items()
    }

    for i, p in enumerate(chunk_list):
        chunk_ds = xr.open_zarr(p).load()

        # sync coord dtypes
        for k in chunk_ds.coords.keys():
            if chunk_ds.coords[k].dtype != coord_dtypes[k]:
                chunk_ds.coords[k] = chunk_ds.coords[k].astype(coord_dtypes[k])

        if i == 0:
            chunk_ds.to_zarr(output_path, mode='w')
        else:
            chunk_ds.to_zarr(output_path, append_dim=region_dim)
    return

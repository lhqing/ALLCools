import numpy as np
import pandas as pd
import pysam
from ALLCools.utilities import reverse_complement, parse_chrom_size
import subprocess


def mode_mc_cov(table, mc, cov):
    mc = table[mc].astype(int)
    cov = table[cov].astype(int)
    frac = np.round(mc / cov).astype(int)
    values = pd.DataFrame({'mc': mc, 'cov': cov, 'frac': frac})
    values = values[['mc', 'cov', 'frac']]
    return values


def mode_mc_uc(table, mc, uc):
    mc = table[mc].astype(int)
    cov = mc + table[uc].astype(int)
    frac = np.round(mc / cov).astype(int)
    values = pd.DataFrame({'mc': mc, 'cov': cov, 'frac': frac})
    values = values[['mc', 'cov', 'frac']]
    return values


def mode_uc_cov(table, uc, cov):
    cov = table[cov].astype(int)
    mc = cov - table[uc].astype(int)
    frac = np.round(mc / cov).astype(int)
    values = pd.DataFrame({'mc': mc, 'cov': cov, 'frac': frac})
    values = values[['mc', 'cov', 'frac']]
    return values


def mode_mc_frac_cov(table, mc_frac, cov):
    frac = table[mc_frac].astype(float)
    cov = table[cov].astype(int)
    mc = np.round(cov * frac).astype(int)
    frac = np.round(frac).astype(int)
    values = pd.DataFrame({'mc': mc, 'cov': cov, 'frac': frac})
    values = values[['mc', 'cov', 'frac']]
    return values


def mode_mc_frac_mc(table, mc_frac, mc):
    frac = table[mc_frac].astype(float)
    mc = table[mc].astype(int)
    cov = np.round(mc / frac).astype(int)
    frac = np.round(frac).astype(int)
    values = pd.DataFrame({'mc': mc, 'cov': cov, 'frac': frac})
    values = values[['mc', 'cov', 'frac']]
    return values


def mode_mc_frac_uc(table, mc_frac, uc):
    frac = table[mc_frac].astype(float)
    uc = table[uc].astype(int)
    cov = np.round(uc / (1 - frac)).astype(int)
    mc = cov - uc
    frac = np.round(frac).astype(int)
    values = pd.DataFrame({'mc': mc, 'cov': cov, 'frac': frac})
    values = values[['mc', 'cov', 'frac']]
    return values


def mode_mc_frac_pseudo_count(table, mc_frac, pseudo_count):
    frac = table[mc_frac].astype(float)
    mc = np.round(pseudo_count * frac).astype(int)
    frac = np.round(frac).astype(int)
    values = pd.DataFrame({'mc': mc, 'frac': frac})
    values['cov'] = pseudo_count
    values = values[['mc', 'cov', 'frac']]
    return values


def get_strand(chrom, pos, fasta):
    base = fasta.fetch(chrom, pos - 1, pos).lower()
    return '+' if base == 'c' else '-'


def get_context(chrom, pos, strand, fasta, num_upstream_bases,
                num_downstream_bases):
    if strand == '+':
        context = fasta.fetch(chrom, pos - 1 - num_upstream_bases,
                              pos + num_downstream_bases).upper()
    else:
        context = reverse_complement(
            fasta.fetch(chrom, pos - 1 - num_downstream_bases,
                        pos + num_upstream_bases).upper())
    return context


def get_strand_and_context(table, chrom, pos, strand, context, fasta_path,
                           num_upstream_bases, num_downstream_bases):
    if strand is not None and context is not None:
        return table.iloc[:, [chrom, pos, strand, context]]

    with pysam.FastaFile(fasta_path) as fasta:
        if strand is None:
            strand = table.shape[1]
            table[strand] = table.apply(
                lambda row: get_strand(row[chrom], row[pos], fasta), axis=1)

        if context is None:
            context = table.shape[1]
            table[context] = table.apply(
                lambda row: get_context(row[chrom], row[pos], row[
                    strand], fasta, num_upstream_bases, num_downstream_bases),
                axis=1)
        return table.iloc[:, [chrom, pos, strand, context]]


def dataframe_to_allc(table,
                      add_chr=False,
                      chrom=0,
                      pos=1,
                      strand=None,
                      context=None,
                      fasta_path=None,
                      chrom_sizes=None,
                      mc=None,
                      uc=None,
                      cov=None,
                      mc_frac=None,
                      pseudo_count=1,
                      num_upstream_bases=0,
                      num_downstream_bases=2):
    # change columns to int
    table.columns = list(range(table.shape[1]))

    # check genome coords
    if (chrom is None) or (pos is None):
        raise ValueError('Must provide chrom and pos')
    if (strand is None) or (context is None):
        if fasta_path is None:
            raise ValueError('Must provide fasta_path if strand or '
                             'context is None, need to use the genome '
                             'FASTA to parse strand or context.')
    if chrom_sizes is not None:
        chroms = set(parse_chrom_size(chrom_sizes).keys())
    elif fasta_path is not None:
        with pysam.FastaFile(fasta_path) as fasta:
            chroms = set(fasta.references)
    else:
        chroms = None

    # check allc values
    if (mc is not None) and (cov is not None):
        value_conversion_func = mode_mc_cov
        mode = 'mc+cov'
    elif (mc is not None) and (uc is not None):
        value_conversion_func = mode_mc_uc
        mode = 'mc+uc'
    elif (uc is not None) and (cov is not None):
        value_conversion_func = mode_uc_cov
        mode = 'uc+cov'
    elif (mc_frac is not None) and (cov is not None):
        value_conversion_func = mode_mc_frac_cov
        mode = 'mc_frac+cov'
    elif (mc_frac is not None) and (mc is not None):
        value_conversion_func = mode_mc_frac_mc
        mode = 'mc_frac+mc'
    elif (mc_frac is not None) and (uc is not None):
        value_conversion_func = mode_mc_frac_uc
        mode = 'mc_frac+uc'
    elif (mc_frac is not None) and (pseudo_count is not None):
        value_conversion_func = mode_mc_frac_pseudo_count
        mode = 'mc_frac+pseudo_count'
    else:
        modes = [
            'mc+cov', 'mc+uc', 'uc+cov', 'mc_frac+cov', 'mc_frac+mc',
            'mc_frac+uc', 'mc_frac+pseudo_count'
        ]
        raise ValueError(f'Need to provide one of these combinations: {modes}')
    # print(f'Using mode {mode} to get cytosine base counts.')

    # add chr to chrom names or not, user specify
    if add_chr:
        table[chrom] = 'chr' + table[chrom].astype(str)
    # select chroms
    if chroms is not None:
        table = table[table[chrom].isin(chroms)].copy()

    # get chrom, pos, strand, context, if needed, parse from fasta
    genome_coords = get_strand_and_context(
        table=table,
        chrom=chrom,
        pos=pos,
        strand=strand,
        context=context,
        fasta_path=fasta_path,
        num_upstream_bases=num_upstream_bases,
        num_downstream_bases=num_downstream_bases)

    # calculate mc, cov, frac (last col)
    if mode == 'mc+cov':
        values = value_conversion_func(table, mc, cov)
    elif mode == 'mc+uc':
        values = value_conversion_func(table, mc, uc)
    elif mode == 'uc+cov':
        values = value_conversion_func(table, uc, cov)
    elif mode == 'mc_frac+cov':
        values = value_conversion_func(table, mc_frac, cov)
    elif mode == 'mc_frac+mc':
        values = value_conversion_func(table, mc_frac, mc)
    elif mode == 'mc_frac+uc':
        values = value_conversion_func(table, mc_frac, uc)
    elif mode == 'mc_frac+pseudo_count':
        values = value_conversion_func(table, mc_frac, pseudo_count)
    else:
        # should have deal with error above
        raise ValueError

    # final allc table
    allc = pd.concat([genome_coords, values], axis=1)
    allc.columns = ['chrom', 'pos', 'strand', 'context', 'mc', 'cov', 'p']
    return allc


def table_to_allc(
        input_path,
        output_prefix,
        sep='\t',
        header=None,
        chunk_size=100000,
        chrom=0,
        pos=1,
        strand=None,
        context=None,
        mc=None,
        uc=None,
        cov=None,
        mc_frac=None,
        pseudo_count=1,
        fasta_path=None,
        num_upstream_bases=0,
        num_downstream_bases=2,
        add_chr=False,
        sort=True
):
    unzip_path = f'{output_prefix}.allc.tsv'
    zip_path = f'{output_prefix}.allc.tsv.gz'
    with open(unzip_path, 'w') as out_f:
        chunks = pd.read_csv(input_path,
                             chunksize=chunk_size,
                             sep=sep,
                             header=header)
        for chunk in chunks:
            # convert table chunks to allc
            allc_chunk = dataframe_to_allc(
                table=chunk,
                add_chr=add_chr,
                chrom=chrom,
                pos=pos,
                strand=strand,
                context=context,
                fasta_path=fasta_path,
                mc=mc,
                uc=uc,
                cov=cov,
                mc_frac=mc_frac,
                pseudo_count=pseudo_count,
                num_upstream_bases=num_upstream_bases,
                num_downstream_bases=num_downstream_bases)
            out_f.write(allc_chunk.to_csv(sep='\t', index=False, header=False))

    # sort, bgzip, tabix allc
    try:
        if sort:
            subprocess.run(
                f'sort -k1,1 -k2,2n {unzip_path} -S 10G | '  # sort
                f'bgzip -c > {zip_path} && '  # bgzip
                f'tabix -f -b 2 -e 2 -s 1 {zip_path} && '  # tabix
                f'rm -f {unzip_path}',  # remove uncompressed file
                shell=True,
                check=True,
                encoding='utf8',
                stderr=subprocess.PIPE)
        else:
            subprocess.run(
                f'bgzip -f {unzip_path} && '  # bgzip
                f'tabix -f -b 2 -e 2 -s 1 {zip_path}',  # tabix
                shell=True,
                check=True,
                encoding='utf8',
                stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(e.stderr)
        raise e
    return zip_path

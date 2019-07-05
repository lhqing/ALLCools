PRESERVED_DIMS = [
    'mc_type',
    'count_type',
    'gene_count_type',
    'strand_type'
]

COORD_ALLOW_VALUES = {
    'count_type': ['mc', 'cov'],
    'gene_count_type': ['exon', 'intron', 'gene_body'],
    'strand_type': ['+', '-', 'both']
}


def check_custom_dim_name(dim_name):
    if dim_name in PRESERVED_DIMS:
        raise ValueError(f'Dimension name {dim_name} is preserved, try another name.')
    if dim_name.startswith('chrom'):
        raise ValueError(f'Dimension name start with "chrom" is preserved, try another name.')
    return

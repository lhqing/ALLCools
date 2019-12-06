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


def check_custom_dim_name_and_return(dim_name):
    if dim_name in PRESERVED_DIMS:
        raise ValueError(f'Dimension name {dim_name} is preserved, try another name.')
    if dim_name.startswith('chrom'):
        raise ValueError(f'Dimension name start with "chrom" is preserved, try another name.')
    if '-' in dim_name:
        _dim_name = dim_name
        dim_name = dim_name.replace('-', '')
        print(f'Dim name {_dim_name} contain "-", will be change to {dim_name} to prevent conflicts. '
              f'Try to provide a name without "-".')
    return dim_name

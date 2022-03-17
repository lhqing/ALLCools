import pandas as pd


def read_gtf(gtf_path):
    gtf = pd.read_csv(gtf_path,
                      comment='#',
                      sep='\t',
                      header=None,
                      names=[
                          'chrom', 'source', 'feature', 'start', 'end',
                          'score', 'strand', 'phase', 'annotation'
                      ])
    return gtf


def subset_gtf(gtf, regions, output_path=None, select_feature=None):
    if (len(regions) == 3) and isinstance(regions[1], int):
        # assume this is a single region
        regions = [regions]

    use_rows = None
    for region in regions:
        chrom, start, end = region
        judge = (gtf['chrom'] == chrom) & (gtf['start'] < end) & (gtf['end'] > start)
        if use_rows is None:
            use_rows = judge
        else:
            use_rows = use_rows | judge
    gtf_sub = gtf[use_rows]

    if select_feature is not None:
        gtf_sub = gtf_sub[gtf_sub['feature'].isin(select_feature)]

    if output_path is not None:
        gtf_sub.to_csv(output_path, sep='\t', index=None, header=None)

    return gtf_sub

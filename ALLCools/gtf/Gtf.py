from collections import defaultdict
from copy import copy
import pathlib

import pandas as pd
from gffutils import FeatureDB, create_db
from gffutils.exceptions import FeatureNotFoundError

from ALLCools.utilities import parse_chrom_size


def create_gtf_db(gtf_path,
                  disable_infer_genes=True,
                  disable_infer_transcripts=True):
    create_db(
        gtf_path, f'{gtf_path}.db',
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts
    )
    return


class Gtf(FeatureDB):
    def __init__(self, *args, **kwargs):
        super(Gtf, self).__init__(*args, **kwargs)

        self._gene_id_to_name = {g.id: g['gene_name'][0]
                                 for g in self.all_features(featuretype='gene')}
        self.gene_ids = set(self._gene_id_to_name.keys())

        # gene name may not be unique
        self._gene_name_to_id_list = defaultdict(list)
        for gene_id, gene_name in self._gene_id_to_name.items():
            self._gene_name_to_id_list[gene_name].append(gene_id)
        return

    def get_gene_name(self, feature_id):
        return self[feature_id].attributes['gene_name'][0]

    def _select_longest_id(self, id_list):
        feature_and_length = {}
        for _id in id_list:
            feature = self[_id]
            length = feature.end - feature.start
            feature_and_length[feature] = length
        # select longest gene id
        _id = sorted(feature_and_length.items(), key=lambda i: i[1])[-1][0]
        return _id

    def get_gene_id_by_name(self, gene_name, select_longest=True):
        try:
            gene_id_list = self._gene_name_to_id_list[gene_name]
        except KeyError:
            raise FeatureNotFoundError(gene_name)

        if len(gene_id_list) == 1:
            gene_id = gene_id_list[0]
        else:
            if select_longest:
                gene_id = self._select_longest_id(gene_id_list)
            else:
                raise ValueError(f'gene_name {gene_name} corresponding to {len(gene_id_list)} '
                                 f'{len(gene_id_list)} transcripts: {gene_id_list}. But '
                                 f'select_longest=False.')
        return gene_id

    def _convert_to_id(self, name):
        if name in self.gene_ids:
            return name
        else:
            return self.get_gene_id_by_name(name)

    def get_gene_feature(self, gene):
        gene_id = self._convert_to_id(gene)
        feature = self[gene_id]
        return feature

    def get_gene_length(self, gene):
        feature = self.get_gene_feature(gene)
        length = feature.end - feature.start
        return length

    def get_gene_transcripts(self, gene, featuretype='transcript'):
        gene_id = self._convert_to_id(gene)
        transcripts = self.children(gene_id, featuretype=featuretype)
        return transcripts

    def get_gene_exons(self, gene, featuretype='exon'):
        gene_id = self._convert_to_id(gene)
        exons = self.children(gene_id, featuretype=featuretype)
        return exons

    def get_gene_tss(self, gene, transcript_featuretype='transcript'):
        transcripts = self.get_gene_transcripts(gene, featuretype=transcript_featuretype)

        # gtf is 1 based
        tss_list = []
        for t in transcripts:
            tss = copy(t)
            if t.strand == '+':
                tss.end = tss.start
            else:
                tss.start = tss.end
            tss.featuretype = 'TSS'
            tss_list.append(tss)
        return tss_list

    def get_gene_promoter(self, gene, chrom_sizes_path, slop=500, transcript_featuretype='transcript'):
        tss_list = self.get_gene_tss(gene, transcript_featuretype=transcript_featuretype)

        chrom_sizes = parse_chrom_size(chrom_sizes_path)

        # gtf is 1 based
        promoter_list = []
        for tss in tss_list:
            promoter = copy(tss)

            promoter.start = max(promoter.start - slop, 0)
            promoter.end = max(promoter.end + slop, chrom_sizes[promoter.chrom])

            promoter.featuretype = 'promoter'
            promoter_list.append(promoter)
        return promoter_list

    def subset_db(self,
                  output_path,
                  genes=None,
                  region_bed=None,
                  span=500000,
                  disable_infer_genes=True,
                  disable_infer_transcripts=True):
        all_features = {}

        if region_bed is not None:
            if isinstance(region_bed, (str, pathlib.Path)):
                region_bed = pd.read_csv(region_bed, sep='\t', index_col=0, header=None)

            for _, (chrom, start, end) in region_bed.iterrows:
                start = max(start - span, 0)
                end = end + span  # it's OK to exceed the chromosome, as here we just select valid gtf features
                features = self.region((chrom, start, end))
                for f in features:
                    all_features[f.id] = f

        if genes is not None:
            for gene in genes:
                feature = self.get_gene_feature(gene)
                start = max(feature.start - span, 0)
                end = feature.end + span
                features = self.region((feature.chrom, start, end))
                for f in features:
                    all_features[f.id] = f

        if len(all_features) == 0:
            print('No features selected.')
            return
        else:
            create_db(list(all_features.values()),
                      dbfn=output_path,
                      disable_infer_genes=disable_infer_genes,
                      disable_infer_transcripts=disable_infer_transcripts)
        return

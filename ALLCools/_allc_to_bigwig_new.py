import pyBigWig
import pysam
from collections import defaultdict
from ._doc import *
from .utilities import parse_mc_pattern, parse_chrom_size


class ContextCounter:
    def __init__(self, mc_contexts):
        self.mc_counts = defaultdict(
            float)  # key is context pattern like CHN, CGN
        self.cov_counts = defaultdict(
            float)  # key is context pattern like CHN, CGN
        # key is each single context like CCC, CGC, value is a list of context pattern this context belong to.
        # A single context may belong to multiple patterns
        self.context_map = defaultdict(list)
        for context in mc_contexts:
            parsed_context_set = parse_mc_pattern(context)
            for c in parsed_context_set:
                self.context_map[c].append(context)

    def add(self, context, mc, cov):
        for c in self.context_map[context]:
            # print(context, c, mc, cov)
            self.mc_counts[c] += mc
            self.cov_counts[c] += cov


class StrandContextCounter:
    def __init__(self, mc_contexts):
        self.mc_watson_counts = defaultdict(
            float)  # key is context pattern like CHN, CGN
        self.cov_watson_counts = defaultdict(
            float)  # key is context pattern like CHN, CGN
        self.mc_crick_counts = defaultdict(
            float)  # key is context pattern like CHN, CGN
        self.cov_crick_counts = defaultdict(
            float)  # key is context pattern like CHN, CGN
        # key is each single context like CCC, CGC, value is a list of context pattern this context belong to.
        # A single context may belong to multiple patterns
        self.context_map = defaultdict(list)
        for context in mc_contexts:
            parsed_context_set = parse_mc_pattern(context)
            for c in parsed_context_set:
                self.context_map[c].append(context)

    def add(self, context, strand, mc, cov):
        # print(context, strand, mc, cov)
        for c in self.context_map[context]:
            if strand == '+':
                self.mc_watson_counts[c] += mc
                self.cov_watson_counts[c] += cov
            else:
                self.mc_crick_counts[c] += mc
                self.cov_crick_counts[c] += cov


def write_entry(counter, context_handle, mc_contexts, strandness, chrom,
                bin_start, bin_size):
    if strandness == 'split':
        for context in mc_contexts:
            # Watson strand
            mc = counter.mc_watson_counts[context]
            cov = counter.cov_watson_counts[context]
            if cov != 0:
                frac = mc / cov
                frac_handle = context_handle[context, '+', 'frac']
                frac_handle.addEntries(chrom,
                                       starts=[bin_start],
                                       values=[frac],
                                       span=bin_size)
                cov_handle = context_handle[context, '+', 'cov']
                cov_handle.addEntries(chrom,
                                      starts=[bin_start],
                                      values=[cov],
                                      span=bin_size)

            # Crick strand
            mc = counter.mc_crick_counts[context]
            cov = counter.cov_crick_counts[context]
            if cov != 0:
                frac = mc / cov
                frac_handle = context_handle[context, '-', 'frac']
                frac_handle.addEntries(chrom,
                                       starts=[bin_start],
                                       values=[frac],
                                       span=bin_size)
                cov_handle = context_handle[context, '-', 'cov']
                cov_handle.addEntries(chrom,
                                      starts=[bin_start],
                                      values=[cov],
                                      span=bin_size)
    else:
        for context in mc_contexts:
            # Both strand
            mc = counter.mc_counts[context]
            cov = counter.cov_counts[context]
            if cov != 0:
                frac = mc / cov
                frac_handle = context_handle[context, 'frac']
                frac_handle.addEntries(chrom,
                                       starts=[bin_start],
                                       values=[frac],
                                       span=bin_size)
                cov_handle = context_handle[context, 'cov']
                cov_handle.addEntries(chrom,
                                      starts=[bin_start],
                                      values=[cov],
                                      span=bin_size)
    return


@doc_params(allc_path_doc=allc_path_doc,
            mc_contexts_doc=mc_contexts_doc,
            chrom_size_path_doc=chrom_size_path_doc,
            strandness_doc=strandness_doc,
            bw_bin_sizes_doc=bw_bin_sizes_doc)
def allc_to_bigwig(allc_path,
                   output_prefix,
                   bin_size,
                   mc_contexts,
                   chrom_size_path,
                   strandness):
    """\
    Generate BigWig files from one ALLC file.

    Parameters
    ----------
    allc_path
        {allc_path_doc}
    output_prefix
        Path prefix of the output BigWig file.
    bin_size
        {bw_bin_sizes_doc}
    mc_contexts
        {mc_contexts_doc}
    strandness
        {strandness_doc}
    chrom_size_path
        {chrom_size_path_doc}
        If chrom_size_path provided, will use it to extract ALLC with chrom order,
        but if region provided, will ignore this.
    """
    if strandness not in {'split', 'both'}:
        raise ValueError(f'strandness need to be "split" or "both", got "{strandness}"')

    chrom_sizes = parse_chrom_size(chrom_size_path)
    chrom_sizes_list = [(k, v) for k, v in chrom_sizes.items()]

    # create bigwig file handles for each case
    # context_handle: key is mC context pattern like CHN, CAN, CGN, value is the output handle
    context_handle = {}
    output_path_collect = {}
    for bw_type in ['frac', 'cov']:
        out_suffix = f'{bw_type}.bw'
        for mc_context in mc_contexts:
            if strandness == 'split':
                file_path = output_prefix + f'.{mc_context}-Watson.{out_suffix}'
                output_path_collect[(mc_context, 'Watson', out_suffix)] = file_path
                # handle for Watson/+ strand
                w_handle = pyBigWig.open(file_path, 'w')
                w_handle.addHeader(chrom_sizes_list)
                context_handle[(mc_context, '+', bw_type)] = w_handle

                file_path = output_prefix + f'.{mc_context}-Crick.{out_suffix}'
                output_path_collect[(mc_context, 'Crick', out_suffix)] = file_path
                # handle for Crick/- strand
                c_handle = pyBigWig.open(file_path, 'w')
                c_handle.addHeader(chrom_sizes_list)
                context_handle[(mc_context, '-', bw_type)] = c_handle
            else:
                # handle for both strand
                file_path = output_prefix + f'.{mc_context}-{strandness}.{out_suffix}'
                output_path_collect[(mc_context, strandness,
                                     out_suffix)] = file_path
                _handle = pyBigWig.open(file_path, 'w')
                _handle.addHeader(chrom_sizes_list)
                context_handle[mc_context, bw_type] = _handle

    def _init_counter(_contexts, _strandness):
        if _strandness == 'split':
            # a counter for +/- strand separately
            _counter = StrandContextCounter(_contexts)
        else:
            # a counter for both +/- strands
            _counter = ContextCounter(_contexts)
        return _counter

    with pysam.TabixFile(allc_path) as allc:
        allc_chroms = set(allc.contigs)
        for chrom, chrom_size in chrom_sizes.items():
            if chrom not in allc_chroms:
                continue
            counter = _init_counter(mc_contexts, strandness)
            cur_bin = 0
            for line in allc.fetch(chrom):
                _, pos, strand, context, mc, cov, _ = line.split('\t')
                pos = int(pos)
                mc = float(mc)
                cov = float(cov)
                this_bin = (pos - 1) // bin_size
                if this_bin != cur_bin:
                    # dump cur_bin counts
                    bin_start = int(cur_bin * bin_size)
                    write_entry(counter=counter,
                                context_handle=context_handle,
                                mc_contexts=mc_contexts,
                                strandness=strandness,
                                chrom=chrom,
                                bin_start=bin_start,
                                bin_size=bin_size)
                    # initiate next bin
                    cur_bin = this_bin
                    counter = _init_counter(mc_contexts, strandness)

                # add counts
                if strandness == 'split':
                    counter.add(context, strand, mc, cov)
                else:
                    counter.add(context, mc, cov)

            # final bin of the chrom
            bin_start = int(cur_bin * bin_size)
            write_entry(counter=counter,
                        context_handle=context_handle,
                        mc_contexts=mc_contexts,
                        strandness=strandness,
                        chrom=chrom,
                        bin_start=bin_start,
                        bin_size=bin_size)
            print(chrom, 'finished')

    for handle in context_handle.values():
        handle.close()
    return output_path_collect

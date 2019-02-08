import os
import contextlib
import gzip
from collections import Counter
from itertools import zip_longest
from .seq import parse_fastq

def hamming(s1, s2):
    '''
    Computes Hamming distance between s1 and s2.
    '''
    if len(s1) != len(s2):
        raise ValueError('{s1} and {s2} must be the same length to compute Hamming distance!'.format(s1=s1, s2=s2))
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1, s2))

def iclip(fq):
    return fq[1][3:7]

def illumina(fq):
    *_, bc = fq1[0].split()[1].split(':')
    return bc


def demux(read1, read2, barcodes, bc_fun=illumina, mismatches=0, data_dir=None, progress=1000000):
    
    fq_tpl='@{}\n{}\n+\n{}\n'
    fname_tpl = 'R{}+{}.fastq.gz'
    barcodes = [s.upper() for s in barcodes]
    output_files = {}
    stat = Counter()
    i = 0
    with contextlib.ExitStack() as stack:
        for fq1,fq2 in zip(parse_fastq(read1), parse_fastq(read2)):
            i += 1
            if i % progress == 0:
                print('{} reads processed'.format(i))
                print(stat)
            bc = bc_fun(fq1)
            for b in barcodes:
                if hamming(b, bc)  <= mismatches:
                    if not b in output_files:
                        output_files[b] = (
                            stack.enter_context(gzip.open(os.path.join(data_dir, fname_tpl.format(1, b)), 'wt')),
                            stack.enter_context(gzip.open(os.path.join(data_dir, fname_tpl.format(2, b)), 'wt'))
                        )
                    output_files[b][0].write(fq_tpl.format(*fq1))
                    output_files[b][1].write(fq_tpl.format(*fq2))
                    stat[b] += 1
                    continue
    return stat

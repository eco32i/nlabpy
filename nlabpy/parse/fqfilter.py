import os
import gzip
import contextlib
import numpy as np
import pandas as pd
from plotnine import *
from nlabpy.parse.seq import parse_fastq
from collections import Counter
from pybloom_live import ScalableBloomFilter

def primer_use(fastq, anchor='GTAATACGACTCACTATAGG', lprimer=25, lrandom=12):
    lanchor = len(anchor)
    pcnt = Counter()
    umis = ScalableBloomFilter()

    cnt = 0
    for seqid, seq, qual in parse_fastq(fastq):
        cnt += 1
        try:
            idx = seq.index(anchor)
            start = idx+lanchor+lrandom
            umi = seq[start-lrandom:start]
            primer = seq[start:start+lprimer]
            if not umi in umis:
                pcnt[primer] += 1
                umis.add(umi)
        except ValueError:
            pcnt['NNNNNNNNNNNNNNNNN'] += 1
        if cnt % 1000000 == 0:
            print('\tprocessed {} records ...\n'.format(cnt))
    return pcnt


def hamming(s1, s2):
    '''
    Computes Hamming distance between s1 and s2.
    '''
    if len(s1) != len(s2):
        raise ValueError('{s1} and {s2} must be the same length to compute Hamming distance!'.format(s1=s1, s2=s2))
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1, s2))


def filter_reads(fastq1, fastq2, pdf, mismatch=1, phase=3,
                   anchor='GTAATACGACTCACTATAGG', lprimer=25, lrandom=12):
    '''
    Allow anchor mismatch
    '''
    
    def _hamming(rec, seq):
        return hamming(rec['primer'], seq)
    
    
    def write_fastq(key, fq1, fq2):
        def format_fastq(fq):
            fastq_tpl = '@{id}\n{seq}\n+\n{qual}\n'
            return fastq_tpl.format(id=fq[0], seq=fq[1], qual=fq[2])
        output_files[key][0].write(format_fastq(fq1))
        output_files[key][1].write(format_fastq(fq2))

    
    
    lanchor = len(anchor)
    # min_len = 4 + lanchor + lrandom + lprimer + 13
    min_len = 4 + lanchor + lrandom + lprimer
    
    primers = set(pdf['primer'])
    dst = pd.DataFrame({'primer': pdf['primer'], 'dist': 0})
    stat = Counter()
    output_files = {}
    umis = ScalableBloomFilter()
    
    with contextlib.ExitStack() as stack:
        # T7 (anchor) sequence not found
        output_files['noanchor'] = (
            stack.enter_context(gzip.open(fastq1.replace('R1', 'R1_noanchor'), 'wt')),
            stack.enter_context(gzip.open(fastq2.replace('R2', 'R2_noanchor'), 'wt'))
        )
        # T7 found but primer is unknown
        output_files['unknown'] = (
            stack.enter_context(gzip.open(fastq1.replace('R1', 'R1_unknown'), 'wt')),
            stack.enter_context(gzip.open(fastq2.replace('R2', 'R2_unknown'), 'wt'))
        )
        # All good
        output_files['good'] = (
            stack.enter_context(gzip.open(fastq1.replace('R1', 'R1_good'), 'wt')),
            stack.enter_context(gzip.open(fastq2.replace('R2', 'R2_good'), 'wt'))
        )
        
        # All primers from T7 sequences
        output_files['T7'] = (
            stack.enter_context(gzip.open(fastq1.replace('R1', 'R1_T7'), 'wt')),
            stack.enter_context(gzip.open(fastq2.replace('R2', 'R2_T7'), 'wt'))
        )
    
        for fqrec1,fqrec2 in zip(parse_fastq(fastq1), parse_fastq(fastq2)):
            stat['total'] += 1
            if len(fqrec1[1]) < min_len or len(fqrec2[1]) < min_len:
                stat['too_short'] += 1
                continue
                
            # Get all possible anchor sequences (due to padding)
            # and their hamming distances
            anchors = {fqrec1[1][i:i+lanchor]: (i, hamming(fqrec1[1][i:i+lanchor], anchor)) 
                       for i in range(phase)}
            idx = None
            for _, d in anchors.items():
                if d[1] < mismatch + 1:
                    stat['anchor'] += 1
                    idx = d[0]
                    break
            
            if idx is not None:    
                start = idx + lanchor + lrandom
                umi = fqrec1[1][idx+lanchor:idx+lanchor+lrandom]
                primer = fqrec1[1][start:start+lprimer]
                
                if not umi in umis:
                    stat['T7'] += 1
                    write_fastq('T7',
                                (fqrec1[0],
                                 fqrec1[1][start:],
                                 fqrec1[2][start:]
                                ),
                                (fqrec2[0],
                                 fqrec2[1][:20],
                                 fqrec2[2][:20]
                                 )
                                )
                    umis.add(umi)
                
                if primer in primers:
                    stat['exact'] += 1
                    write_fastq('good', fqrec1, fqrec2)
                else:
                    dst['dist'] = dst.apply(_hamming, axis=1, args=(primer,))
                    if dst.loc[dst['dist'].argmin(), 'dist'] < 2:
                        # Inexact primer matches still go to the "good" file
                        stat['inexact'] += 1
                        write_fastq('good', fqrec1, fqrec2)
                    else:
                        # No primer found
                        stimport os
import gzip
import contextlib
import numpy as np
import pandas as pd
from plotnine import *
from nlabpy.parse.seq import parse_fastq
from collections import Counter
from pybloom_live import ScalableBloomFilterat['unknown'] += 1
                        write_fastq('unknown', fqrec1, fqrec2)    
            else:
                # T7 (anchor) not found
                stat['NA'] += 1
                write_fastq('noanchor', fqrec1, fqrec2)
            if stat['total'] % 1000000 == 0:
                print('\tprocessed {} records ...\n'.format(stat['total']))
                print(stat)
    #print(stat)
    return stat

def read_stat(fastq1, fastq2, pdf, mismatch=1, phase=3,
                   anchor='GTAATACGACTCACTATAGG', lprimer=25, lrandom=12):
    '''
    Allow anchor mismatch
    '''
    
    def _hamming(rec, seq):
        return hamming(rec['primer'], seq)
    
    
    lanchor = len(anchor)
    min_len = 2 + lanchor + lrandom + lprimer
    primers = set(pdf['primer'])
    
    dst = pd.DataFrame({'primer': pdf['primer'], 'dist': 0})
    reads = []
    
    for fqrec1,fqrec2 in zip(parse_fastq(fastq1), parse_fastq(fastq2)):
        read = {'R1_len': len(fqrec1[1]),
               'R2_len': len(fqrec2[1]),
               'padding': -1,}
        if len(fqrec1[1]) < min_len:
            continue
        
        anchors = {fqrec1[1][i:i+lanchor]: (i, hamming(fqrec1[1][i:i+lanchor], anchor)) 
                   for i in range(phase)}
        idx = None
        for _,d in anchors.items():
            if d[1] < mismatch + 1:
                stat['anchor'] += 1
                idx = d[0]
                break
        if idx is not None:
            start = idx + lanchor + lrandom
            primer = fqrec1[1][start:start+lprimer]
            umi = fqrec1[1][idx+lanchor:idx+lanchor+lrandom]
                
            read['padding'] = idx

            if primer in primers:
                read['status'] = 'exact'
                read['umi'] = umi
                read['pstart'] = pdf[pdf['primer']==primer]['position'].values[0]
            else:
                dst['dist'] = dst.apply(_hamming, axis=1, args=(primer,))
                if dst.loc[dst['dist'].argmin(), 'dist'] < 2:
                    read['status'] = 'inexact'
                    read['umi'] = umi
                    read['pstart'] = pdf[pdf['primer']==dst.loc[dst['dist'].argmin(), 'primer']]['position'].values[0]
                else:
                    read['status'] = 'unknown'
                    read['pstart'] = 0
        else:
            read['status'] = 'NA'
            read['pstart'] = 0
        reads.append(read)
        
    return pd.DataFrame.from_records(reads)

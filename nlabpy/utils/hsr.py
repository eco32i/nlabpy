import pandas as pd
import numpy as np
import sys

# The following code is taken from Ben Langmead's lab teaching website:
# https://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_kEditDp.ipynb
# and modified slightly

# Assume x is the string labeling rows of the matrix and y is the
# string labeling the columns

def trace(D, x, y):
    ''' Backtrace edit-distance matrix D for strings x and y '''
    i, j = len(x), len(y)
    xscript = []
    while i > 0:
        diag, vert, horz = sys.maxsize, sys.maxsize, sys.maxsize
        delt = None
        if i > 0 and j > 0:
            delt = 0 if x[i-1] == y[j-1] else 1
            diag = D[i-1, j-1] + delt
        if i > 0:
            vert = D[i-1, j] + 1
        if j > 0:
            horz = D[i, j-1] + 1
        if diag <= vert and diag <= horz:
            # diagonal was best
            xscript.append('R' if delt == 1 else 'M')
            i -= 1; j -= 1
        elif vert <= horz:
            # vertical was best; this is an insertion in x w/r/t y
            xscript.append('I')
            i -= 1
        else:
            # horizontal was best
            xscript.append('D')
            j -= 1
    # j = offset of the first (leftmost) character of t involved in the
    # alignment
    return j, (''.join(xscript))[::-1] # reverse and string-ize

def kEditDp(p, t):
    ''' Find and return the alignment of p to a substring of t with the
        fewest edits.  We return the edit distance, the offset of the
        substring aligned to, and the edit transcript.  If multiple
        alignments tie for best, we report the leftmost. '''
    D = np.zeros((len(p)+1, len(t)+1), dtype=int)
    # Note: First row gets zeros.  First column initialized as usual.
    D[1:, 0] = range(1, len(p)+1)
    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            delt = 1 if p[i-1] != t[j-1] else 0
            D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
    # Find minimum edit distance in last row
    # Replace the original loop with numpy functions
    #
    # mnJ, mn = None, len(p) + len(t)
    #for j in range(len(t)+1):
    #    if D[len(p), j] < mn:
    #        mnJ, mn = j, D[len(p), j]
    mnJ = np.argmin(D[len(p),:])
    mn = np.min(D[len(p),:])
    # Backtrace; note: stops as soon as it gets to first row
    off, xcript = trace(D, p, t[:mnJ])
    # Return edit distance, offset into T, edit transcript
    return mn, off, xcript, D

def rc(s):
    '''
    Reverse complement to DNA sequence given in s
    '''
    NUC = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    _rc = [NUC[x] for x in s]
    return ''.join(_rc)[::-1]

def hamming(s1, s2):
    '''
    Computes Hamming distance between s1 and s2.
    '''
    if len(s1) != len(s2):
        raise ValueError('{s1} and {s2} must be the same length to compute Hamming distance!'.format(s1=s1, s2=s2))
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1, s2))


def cmp_reads(fastq, pdf, ref, queue=None):
    '''
    Allow anchor mismatch
    '''
    
    result = []
    cnt = 0
    for fqrec in parse_fastq(fastq):
        start = 20 + 12 # lanchor + lrandom
        #umi = fqrec1[1][idx+lanchor:idx+lanchor+lrandom]
        primer = fqrec[1][start:start+25]
        read = fqrec[1][start+25:]
        
        if pdf[pdf['primer']==primer]['position'].values:
            pstart = pdf[pdf['primer']==primer]['position'].values[0]
            ref_read = ref[pstart-len(read):pstart]
        else:
            cnt += 1
            continue
        if len(read) == len(ref_read):
            dist = hamming(rc(read), ref_read)
        else:
            dist = -1
        result.append({
            'pstart': pstart,
            'read': read,
            'rlen': len(read),
            'ref': ref_read,
            'dist': dist
        })
    print(cnt)
    return pd.DataFrame.from_records(result)

def cmp_reads_DP(fastq, pdf, ref, queue=None):
    '''
    Allow anchor mismatch
    '''
    
    result = []
    cnt = 0
    for fqrec in parse_fastq(fastq):
        start = 20 + 12 # lanchor + lrandom
        #umi = fqrec1[1][idx+lanchor:idx+lanchor+lrandom]
        primer = fqrec[1][start:start+25]
        read = fqrec[1][start+25:]
        
        if pdf[pdf['primer']==primer]['position'].values:
            pstart = pdf[pdf['primer']==primer]['position'].values[0]
            ref_read = ref[pstart-len(read)-1:pstart]
        else:
            cnt += 1
            continue
        dist, off, xscript, *_ = kEditDp(rc(read), ref_read)
        
        result.append({
            'pstart': pstart,
            'read': read,
            'rlen': len(read),
            'ref': rc(ref_read),
            'dist': dist,
            'offset': off,
            'xs': xscript
        })
    print(cnt)
    # return cnt of reads with no primer and df for the rest of reads
    return cnt, pd.DataFrame.from_records(result)

def mutations_per_position(fastq, ref, pdf, start=25, mismatch=4):
    ref_len = len(ref)
    df = pd.DataFrame({
        'ref': list(ref),
        'pos': range(1,ref_len+1),
        'A': ref_len*[0],
        'G': ref_len*[0],
        'C': ref_len*[0],
        'T': ref_len*[0],
        'I': ref_len*[0],
        'D': ref_len*[0],
        'N': ref_len*[0],
        })
    for fqrec in parse_fastq(fastq):
        primer = fqrec[1][start:start+25]
        read = fqrec[1][start+25:]

        if pdf[pdf['primer']==primer]['position'].values:
            pstart = pdf[pdf['primer']==primer]['position'].values[0]
            ref_read = ref[pstart-len(read)-1:pstart]
        else:
            continue
        dist, off, xscript, *_ = kEditDp(rc(read), ref_read)
        if dist < mismatch:
            for i, pos in enumerate(range(pstart-len(read)-1, pstart-1)):
                df.at[pos, rc(read)[i]] += 1
    df['counts'] = df[['A', 'G', 'C', 'T']].sum(axis=1)
    return df


def rlength_per_position(fastq, ref, pdf, start=25, mismatch=4):
    res = []

    for fqrec in parse_fastq(fastq):
        primer = fqrec[1][start:start+25]
        read = fqrec[1][start+25:]

        if pdf[pdf['primer']==primer]['position'].values:
            pstart = pdf[pdf['primer']==primer]['position'].values[0]
            ref_read = ref[pstart-len(read)-1:pstart]
        else:
            continue
        dist, off, xscript, *_ = kEditDp(rc(read), ref_read)
        if dist < mismatch:
            res.append({'pos': pstart, 'length': len(read)})

    return pd.DataFrame.from_records(res)

from bokeh.models import ColumnDataSource
import pandas as pd
import pysam

def chrom_stat(bamfile, chromosomes=None):
    '''
    Given a bamfile, counts the number of reads in each chromosome
    (or rather a reference sequence as determined by bam header)
    
    Returns a dataframe with chromosome names, sizes, and midpoint coords.
    '''
    
    def count_chr_reads(rec):
        cnt = 0
        chr_iter = bam.fetch(rec['SN'], 0, rec['LN'])
        for _ in chr_iter:
            cnt += 1
        return cnt
    
    bam = pysam.AlignmentFile(bamfile, 'rb')
    df = pd.DataFrame.from_records(bam.header['SQ'])
    if chromosomes:
        df = df[df['SN'].isin(chromosomes)].copy()
    df['reads'] = df.apply(count_chr_reads, axis=1)
    df['center'] = df['LN'] // 2
    return df


def bam_to_chromosomes(bamfile, bins=100, colorspec='rgba(255, 0, 0, {:.2f})', chromosomes=None, reads=None):
    '''
    Divides chromosomes into specified number of bins and counts number of
    reads in each bin.
    
    Sets the color (given in a colorspec string) opacity proportional to the number of reads.
    '''
    
    def get_color_spec(rec):
        value = rec['read_count'] / reads_max
        return colorspec.format(value)

    bam = pysam.AlignmentFile(bamfile, 'rb')
    df_chr = chrom_stat(bamfile, chromosomes)
    if chromosomes is None:
        chromosomes = df_chr['SN']
    res = []
    for chro in chromosomes:
        bin_width = int(df_chr[df_chr['SN'] == chro]['LN'] / bins)
        for i in range(bins):
            cnt = 0
            start = bin_width * i
            end = start + bin_width
            chriter = bam.fetch(chro, start, end)
            for r in chriter:
                if reads is not None:
                    if r.qname in reads:
                        cnt += 1
                else:
                    cnt += 1
            res.append((start, end, cnt, chro ))

    df = pd.DataFrame.from_records(res, columns=['start', 'end', 'read_count', 'chro'])
    df['width'] = df.end - df.start
    df['center'] = df.start + df.width // 2
    
    reads_max = df.read_count.max()
    df['colorspec'] = df.apply(get_color_spec, axis=1)
    return df


def chrom_to_int(data, values, field="chro"):
    df = data.copy()
    df[field] = 0
    for i,val in enumerate(reversed(values)):
        df.ix[df['SN']==val, field] = i
    return df
        

def draw_chromosomes(p, data, chromosomes=None):
    df = chrom_to_int(data, chromosomes)
    p.rect(x="center", y="chro", width="LN", height=0.6, source=ColumnDataSource(df),
          fill_alpha=0, line_alpha=0.4)
    
    
def draw_ideograms(p, data, chromosomes=None):
    df = chrom_to_int(data, chromosomes)

    for stain,group in df.groupby('gieStain'):
        if stain == 'acen':
            starts = group.start.values.tolist()
            ends = group.end.values.tolist()
            chros = group.chro.values.tolist()
            xcoords_left = [
                (start, end, start) for start,end in 
                zip(starts[0::2], ends[0::2])
            ]
            xcoords_right = [
                (start, end, end) for start,end in 
                zip(starts[1::2], ends[1::2])
            ]
            ycoords_left = [(chro+0.3, chro, chro-0.3) for chro in chros[0::2]]
            ycoords_right = [(chro, chro+0.3, chro-0.3) for chro in chros[1::2]]
            xs = [x for coords in zip(xcoords_left, xcoords_right) for x in coords]
            ys = [x for coords in zip(ycoords_left, ycoords_right) for x in coords]
            p.patches(
                xs=xs,
                ys=ys,
                fill_alpha=0.5, line_alpha=0.1, color="#2244FF")
        else:
            src = ColumnDataSource(group)
            p.rect(x="center", y="chro", width="width", height=0.6, source=src,
                   fill_alpha=0.15, line_alpha=0, color="color")
    
    
def draw_track(p, data, chromosomes=None):
    df = chrom_to_int(data, chromosomes)
    p.rect(x="center", y="chro", width="width", height=0.8, source=ColumnDataSource(df),
           fill_alpha=1, line_alpha=0, color="colorspec")    


def draw_names(p, chromosomes):
    p.text(
        [-1500000 for x in chromosomes],
        [i - 0.5 for i,_ in enumerate(reversed(chromosomes))],
        text=list(reversed(chromosomes)),
        text_align="right")

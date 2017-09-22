import numpy as np
from bokeh.plotting import vplot
from bokeh.palettes import brewer
from bokeh.models import FixedTicker, NumeralTickFormatter, Range1d, LinearAxis


def strand_data(bamfiles, reference, region):
    '''
    Returns a dataframe with coverage data by strand for the given bamfiles,
    reference and region
    
    `bamfiles`: dictionary of the form {'plus': <filename>, 'minus': <filename>}
        where filename refers to the corresponding bamfile (strand-split)
    
    `reference`: reference, as specified in the bamfile header
    
    `region`: tuple (min, max) of the positions that define the region
    '''
    cov = []
    for strand,filename in bamfiles.items():
        bamfile = pysam.AlignmentFile(filename, 'rb')
        start,end = region
        pileup = bamfile.pileup(reference, start, end,
                               truncate=True,
                               max_depth=1e9)
        for col in pileup:
            cov.append((col.reference_pos, col.n, strand))
    df = pd.DataFrame.from_records(cov, columns=('position', 'coverage', 'strand'))
    df['coverage'] = np.where(df['strand'] == 'minus', -1*df['coverage'], df['coverage'])
    return df



def strand_coverage(strands, limits, title, anno=None):
    _strands = ['plus', 'minus'] 
    data = strand_data(strands, reference, limits)
    xs = []
    ys = []
    for strand in _strands:
        region = list(data[data['strand'] == strand]['position'])
        xs.append(np.hstack((region, region[::-1])))
        ys.append(np.hstack((data[data['strand']==strand]['coverage'], np.zeros(len(region)))))

    colors = brewer["Spectral"][4]
    alphas = [0.8, 0.5]
    max_yval = max([np.max(np.abs(y)) for y in ys])
    max_xval = max([np.max(x) + 100 for x in xs])
    min_xval = min([np.min(x) for x in xs])
    p = figure(title=title,
               y_range=(-max_yval, max_yval),
               x_range=(min_xval, max_xval+1))

    p.patches(xs, ys, color="#268bd2",
             alpha=alphas, line_color=None)
    
    if anno is not None:
        tracks = set(anno['track'])
        anno['trk'] = len(tracks) + 0.5
        anno.ix[anno.strand=='+', 'trk'] *= -1
        anno['trkname'] = anno['trk']
        for i,track in enumerate(tracks):
            anno.ix[(anno.track==track) & (anno.strand=='-'), 'trk'] -= i / 2
            anno.ix[(anno.track==track) & (anno.strand=='-'), 'trkname'] -= (i/2) + 0.25
            anno.ix[(anno.track==track) & (anno.strand=='+'), 'trk'] += i / 2
            anno.ix[(anno.track==track) & (anno.strand=='+'), 'trkname'] += (i/2) - 0.25
        p.extra_y_ranges = {'tracks': Range1d(start=-len(tracks)-1, end=len(tracks)+1)}
        p.add_layout(LinearAxis(y_range_name='tracks'), 'right')
        
        p.rect(x='center', y='trk', width='length', height=0.25, fill_color='blue',
              source=ColumnDataSource(anno), fill_alpha=0.2, y_range_name='tracks')
        p.text(x='center', y='trkname', text='name', text_align='center',
               source=ColumnDataSource(anno), y_range_name='tracks',
              text_font_size='10pt', text_baseline="middle")
        p.text(x=max_xval, y='trk', text='track', text_align='right',
               source=ColumnDataSource(anno), y_range_name='tracks', text_baseline='middle')
        
        mark_dist = (max_xval - min_xval) // 20
        marks = [[(i, ann['trk'], ann['strand']) for i in range(min_xval, max_xval, mark_dist)
                if ann['start'] < i < ann['end']] for _,ann in anno.iterrows()]
        
        for m in marks:
            p.text(x=[x[0] for x in m],
                   y=[x[1] for x in m],
                   text=['>' if x[2] == '+' else '<' for x in m],
                   text_align='right',
                   y_range_name='tracks', text_baseline='middle')
            

    p.plot_width = 900
    p.grid.grid_line_color = None
    p.yaxis.minor_tick_line_color = None
    p.yaxis[1].major_tick_line_color = None
    p.yaxis[1].major_label_text_color = None
    p.xaxis.minor_tick_line_color = None
    p.xaxis[0].ticker = FixedTicker(ticks=[limits[0], limits[0]+(limits[1]-limits[0])//2, limits[1]])
    p.xaxis[0].formatter = NumeralTickFormatter(format="0,0")
    p.yaxis[0].formatter = NumeralTickFormatter(format="(0,0)")
    return p


def region_plot(feature, padding=50, limits=None):
    
    def _fcheck(rec):
        s = rec['TU'].replace('_', '-')
        return feature in s.split('-')
    
    vg = gff[gff['gene'] == feature][['start','end']].values
    vu = utr[utr.apply(_fcheck, axis=1)][['start','end']].values
    
    if limits is None:
        region = np.vstack((vg, vu))
        rmax = np.max(region)
        rmin = np.min(region)
        if padding is not None:
            rmax += padding
            rmin -= padding
        limits = (rmin, rmax)
    
    utrs = utr[utr.apply(_fcheck, axis=1)].copy()
    utrs['length'] = utrs['end'] - utrs['start']
    utrs.sort_values('length', ascending=False, inplace=True)
    
    utr_= utrs.iloc[0].to_frame().transpose()
    utr_['name'] = feature
    utr_['track'] = 'UTR'
    gene = gff[gff['gene'] == feature].copy()
    gene['track'] = 'gene'
    gene['name'] = feature
    
    genes = pd.concat([gene, utr_])
    genes['center'] = genes.start + (genes.end - genes.start) // 2
    genes['length'] = genes.end - genes.start

    p1 = strand_coverage(strands_nobcm, limits, "WT -BCM", genes)
    p2 = strand_coverage(strands_bcm, limits, "WT +BCM", genes)
    p = vplot(p1, p2)
    show(p)

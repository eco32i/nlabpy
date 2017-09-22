import pandas as pd
import numpy as np

def extract_gff(line, ftype='gene'):
    rec = {}
    fields = line.strip().split('\t')
    if fields[2] != ftype:
        return None
    opts = {}
    for f in fields[8].split(';'):
        k,v = f.split('=')
        opts[k] = v
        
    rec.update({
            'gene': opts['gene'],
            'start': int(fields[3]),
            'end': int(fields[4]),
            'strand': fields[6]
        })
    return rec


def extract_bed(line):
    rec = {}
    fields = line.strip().split('\t')
    rec.update({
            'TU': fields[3],
            'start': int(fields[1]),
            'end': int(fields[2])
        })
    return rec



def anno_df(annofile, extract=extract_gff):
    result = []
    with open(annofile) as fi:
        for line in fi:
            if not line.startswith('#'):
                rec = extract(line)
                if rec:
                    result.append(rec)
    return pd.DataFrame.from_records(result)


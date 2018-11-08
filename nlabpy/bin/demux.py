#!/usr/bin/env python
import argparse
from nlabpy.parse.demux import demux

def main():
    parser = argparse.ArgumentParser(
        description='demux: demultiplexing Illumina .fastq files.')

    parser.add_argument('R1', type=str, help='.fastq file with read1 sequences')
    parser.add_argument('R2', type=str, help='.fastq file with read2 sequences')
    parser.add_argument('-b', '--barcodes', nargs='+', type=str, help='barcode sequences')
    parser.add_argument('-m', '--mismatches', type=int, default=0, help='maximum number of mismatches in barcode')
    parser.add_argument('-d', '--data_dir', type=str, default='',
            help='directory to write output files to')

    args = parser.parse_args()
    kwargs = vars(args)
    
    print(args.R1)
    print(args.R2)
    print(kwargs)

    stat = demux(args.R1, args.R2, args.barcodes,
            mismatches=kwargs['mismatches'], data_dir=kwargs['data_dir'])
    print(stat)
    print('Done.')

# Run the actual script here.
    
if __name__ == "__main__":
    main()

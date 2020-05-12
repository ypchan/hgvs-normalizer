#!/usr/bin/env python3

'''
This script is written to parse fasta file into json file.
Author: chenyanpeng@dnastories
Date:   2020-05-08
'''

import sys
import json
import gzip
import argparse

def fasta2dict(fasta_file):
    fa_dict = dict()

    fa_fh = gzip.open(fasta_file,'rt') if fasta_file.endswith('.gz') else open(fasta_file,'rt')
    for line in fa_fh:
        line = line.rstrip('\n')
        if line.startswith(">"):
            identifer = line.split()[0].lstrip('>')
            fa_dict[identifer] = []
        else:
            fa_dict[identifer].append(line)
    fa_dict = {identifier:''.join(seq) for identifier, seq in fa_dict.items()}
    fa_fh.close()
    return fa_dict

def parse_args(argv):
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="Example: genome2json.py file.fasta file.json -g")

    parser.add_argument("infile", metavar='file.fasta', type=str,
            help='input file in FASTA format (compressed gz file is allowed)')
    parser.add_argument("outfile", metavar='file.json', type=str,
            help='output file in JSON format')
    parser.add_argument('-g', '--gzip', action='store_true',
            help='gzip compress the out file')

    args = parser.parse_args()
    return args

def outjson(fa_dict, outfile, gztag):
    if gztag:
        outfile = outfile + '.gz'
        out_f = gzip.open(outfile, 'wt')
    else:
        out_f = open(outfile, 'wt')
    json.dump(fa_dict, out_f, indent=2)
    out_f.close()

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    fa_dict = fasta2dict(args.infile)
    outjson(fa_dict, args.outfile, args.gzip)
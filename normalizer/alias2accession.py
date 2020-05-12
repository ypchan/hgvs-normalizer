#!/usr/bin/env python3
'''
This script is written to parse the accession of chr and its alias into json file.
Author : chenyanpeng@dnastories
Date   : 2020-05-08
Version: 1.0

Format: inputfile.tsv
6   NC11112.1
9   NC12344.2
'''

import sys
import json
import fileinput
import argparse

def parse_args(argv):
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            usage='\n  %(prog)s [-h] <id2chr.tsv> <id2chr.json> [-k 1] [-v 2]',
            epilog="example:\n  %(prog)s id2chr.tsv id2chr.json -k 1 -v 2")

    parser.add_argument("infile", metavar='id2chr.tsv', type=str,
            help='input file in TSV format, the first field is the alias of accession (second field) of chr')
    parser.add_argument("outfile", metavar='id2chr.json', type=str,
            help='output file in JSON format')
    parser.add_argument("-k", '--key_field', metavar='int', type=int, default=1,
            help='set the key field number [default: 1]')
    parser.add_argument("-v", '--value_field', metavar='int', type=int, default=2,
            help='set the value field number [default: 2]')
    args = parser.parse_args()
    return args

def get_idxchr_dict(idchr_file, key_field, value_field):
    idxchr_dict = dict()
    with fileinput.input(files=idchr_file) as infh:
        for line in infh:
            line = line.rstrip('\n')
            line_lst = line.split('\t')
            key = line_lst[key_field - 1]
            value = line_lst[value_field - 1]
            idxchr_dict[key] = value
    return idxchr_dict

def outjson(idxchr_dict, outfile):
    with open(outfile, 'wt') as outf:
        json.dump(idxchr_dict, outf, indent=2)

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    idxchr_dict = get_idxchr_dict(args.infile, args.key_field, args.value_field)
    outjson(idxchr_dict, args.outfile)
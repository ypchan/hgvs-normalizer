#!/usr/bin/env python3
'''
This script is written to format input TSV file containing the variants information.
Author : chenyanpeng@dnastories
Date   : 2020-05-11
Version: 1.0
------------------------------------------------------
Input: variants.tsv
6       g.55146556_55146557insC
14      g.13687521_13687522insATT
14      g.13742402A>C
14      g.13584928_13584929insC

Output: variants.tsv.formatted
NC_006588.3:g.55146556_55146557insC
NC_006596.3:g.13687521_13687522insATT
NC_006596.3:g.13742402A>C
NC_006596.3:g.13584928_13584929insC
------------------------------------------------------
'''

import sys
import json
import gzip
import argparse
import fileinput

def parse_args():
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            usage='\n  %(prog)s [-h] <variants.tsv> <id2accession.json>',
            epilog="example:\n  %(prog)s variants.tsv id2accession.json")

    parser.add_argument("tsvfile", metavar='variants.tsv', type=str,
            help='input file in TSV format')
    parser.add_argument("jsonfile", metavar='id2accession.json', type=str,
            help='JSON file recording the accession and its alias id')
    
    args = parser.parse_args()
    return args

def load_json(jsonfile):
	try:
		with open(jsonfile, 'r') as jsonfile_f:
    		json_dict = json.load(jsonfile_f)
    except:
    	sys.exit(f'Error: failed to open {jsonfile}')
    return json_dict

def formatting(variant_record, json_dict):
	variant_record = variant_record.rstrip('\n')
	record_lst = variant_record.split('\t')
	num_chr = record_lst[0]
	accession = json_dict.get(num_chr)
	if not accession:
		sys.exit(f'{num_chr} is not in jsonfile')
	return accession + ':' + record_lst[1]

def output(tsvfile, json_dict):
	with fileinput.input(files=tsvfile) as tsv_fh, open(tsvfile + '.formatted', 'wt') as outfh:
		for line in tsv_fh:
			record_formatted = formatting(line, json_dict)
			outfh.write(record_formatted + '\n')

if __name__ == '__main__':
	args = parse_args()
	json_dict = load_json(args.jsonfile)
	output(args.tsvfile, json_dict)
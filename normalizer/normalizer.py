#!/usr/bin/env python3

'''
This script is written to normalize the mutation description.
Author: chenyanpeng@dnastories
Date  : 2020-05-08
'''

import re
import sys
import json
import argparse

def parse_args(argv):
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="Example: normalizer.py variant.tsv variant.tsv.noramlized variant.tsv.error")

    parser.add_argument("infile", metavar='variant.tsv', type=str,
            help='input file in TSV format')
    parser.add_argument('-o', '--outfile', metavar='<variant.tsv.noramlized>', type=str, default='variant.tsv.noramlized',
            help='output file in TSV format')
    parser.add_argument('-e', '--errorfile', metavar='<variant.tsv.error>', type=str, default='variant.tsv.error',
            help='records that failed to be normalize')
    
    args = parser.parse_args()
    return args

def load_json(genomejson, rnajson, proteinjson):
	# JSON file is similar to dictionary in Pyton
	try:
		with open(genomejson, 'r') as genomejson_f:
    		genome_dict = json.load(genomejson_f)
    	with open(rnajson, 'r') as rnajson_f:
    		rna_dict = json.load(rnajson_f)
    	with open(proteinjson, 'r') as proteinjson_f:
    		protein_dict = json.load(proteinjson_f)
    	return genome_dict, rna_dict, protein_dict
    except IOError:
    	sys.exit('Error: load json file error.')

def get_g_varstring(variants_tsvfile:str, g_varstring_field:int) -> object:
	with open(variants_tsvfile, 'rt') as tsv_file:
		for line in tsv_file:
			line = line.rstrip('\n')
			line_lst = line.split('\t')
			line_lst = [i.strip() for i in line_lst]
			g_varstring = line_lst[int(g_varstring_field) - 1]
			yield g_varstring

def parse_variant(g_varstring:str, genome_json_dict:dict, output_lst:list, error_file_fh:object) -> list:
	acc_level_posedit_lst = re.split(r'[:.]', g_varstring)
	# acc_level_posedit_lst = [accession, level, posedit]

	if len(acc_level_posedit_lst) != 3:
		error_file_fh.write(g_varstring + '\n')
		sys.exit(f'unknown {g_varstring}')
		
	accession = acc_level_posedit_lst[0]
	seq = genome_json_dict.get(accession)
	if not seq:
		error_file_fh.write(g_varstring + '\n')
		sys.exit(f'Error: faild to find {g_varstring} in genome json file')

	posedit = acc_level_posedit_lst[2]
	if '>' in posedit:
		output_list.append(g_varstring)
		return None
	elif 'delins' in posedit:
		output_list.append(g_varstring)
		return None
	elif 'inv' in posedit:
		output_list.append(g_varstring)
		return None
	elif 'con' in posedit:
		output_list.append(g_varstring)
		return None
	elif 'dup' in posedit:
		return ['dup', acc_level_posedit_lst]
	elif 'del' in posedit:
		return ['del', acc_level_posedit_lst]
	elif 'ins' in posedit:
		return ['ins', acc_level_posedit_lst]
	else:
		error_file_fh.write(g_varstring + '\n')
		sys.exit(f'Error: invalid g_varstring {posedit}')

def accession2seq(accession:str, var_level:str, json_lst: lst) -> str:
	genome_dict, rna_dict, protein_dict = json_lst
	if var_level == 'r':
		seq = rna_dict.get(accession)
	elif var_level == 'p':
		seq = protein.get(accession)
	else:
		seq = genome_dict.get(accession)
	return seq

def shift_ins_to_3end(seq:str, posedit:str):
	locus = re.search(r'(\d+)_(\d+)ins', posedit)
	num_left =  int(locus.group(1)) 
	num_right = int(locus.group(2)) 

	ref_from_string = posedit.split('ins')[1]
	if not ref_from_string:
		sys.exit(f'Error: insertion site {posedit} no insertion sequence', file=sys.stderr, flush=True)
		ref_tag = False
	else:
		ref_tag = True

	# before insertion
	seq_left10 = seq[num_left - 10: num_left]
	seq_right10 = seq[num_right - 1: num_right + 9]
	seq20_before_ins = seq_left10 + seq_right10
	# after insertion
	seq20_after_ins = seq_left10 + ref_from_string + seq_right10

	for i in range(11, 20):
		shift_seq20_after_ins = seq_ref20[:i] + ref_from_string + seq_ref20[i:]
		if shift_seq20_after_ins != seq20_after_ins:
			break
		shift = True
	locus = [num_left + i - 11, num_right + i - 11 ]
	if shift:
		print(f'old: {num_left}_{num_right} -> new: {locus[0]}_{locus[1]}', sys.stdout, flush=True)
	
	g_posedit = '_'.join([str(pos) for pos in locus]) + 'ins' + ref_from_string
	return  g_posedit

def shift_del_to_3end(seq:str, posedit:str) -> str:
	if '_' in posedit: # deletion of >=2 bases 
		try:
			locus = re.search(r'(\d+)_(\d+)del', posedit)
			num_left =  int(locus.group(1)) 
			num_right = int(locus.group(2)) 
			num_bases = True
		except:
			sys.exit(f'Error: parse {posedit} error')
	else:              # single base deletion 
		try:
			locus = re.search(r'(\d+)del', posedit)
			num_left = int(locus.group(1))
			num_right = num_left
			num_bases = False
		except:
			sys.exit(f'Error: parse {posedit} error')
	
	ref_from_locus = seq[num_left - 1:num_right] # bases of deletion by indexing 
	ref_from_string = posedit.split('del')[1] # bases of deletion from posedit
	if ref_from_string != ref_from_locus:
		sys.exit(f'Error: {posedit} length not equal to end - start + 1', file=sys.stderr, flush=True)
	
	# before deletion
	seq_left10 = seq[num_left - 11: num_left - 1]
	seq_right10 = seq[num_right : num_rigth + 10]
	deletion = seq[num_left - 1: num_right]
	seq20_before_del = seq_left10 + deletion + seq_right10
	# after deletion
	seq_after_deletion = seq_left10 + seq_right10
	
	step = len(deletion)
	len_seq = step + 20
	for i in range(10 + step, len_seq, step):
		shift_deletion = seq[i: i + step]
		if shift_deletion != deletion:
			break
		shift = True

	locus = [num_left + i - 10, num_left + i - 10 + step - 1]
	if shift:
		if num_bases:
			print(f'old: {num_left}_{num_right} -> new: {locus[0]}_{locus[1]}', file=sys.stdout, flush=True)
		else:
			print(f'old: {num_left} -> new: {locus[0]}', file=sys.stdout, flush=True)

	if num_bases:
		posedit = '_'.join([str(pos) for pos in locus]) + 'del' + ref
	else:
		posedit = ''.join([str(locus[0]), 'del', ref])
	return g_posedit 

def shift_dup_to_3end(seq:str, posedit:str) -> str:
	if '_' in posedit: # duplication of >= 2 bases
		try:
			locus = re.search(r'(\d+)_(\d+)dup', posedit)
			num_left =  int(locus.group(1)) 
			num_right = int(locus.group(2)) 
			num_bases = True
		except:
			sys.exit(f'Error: parse {posedit} error')
	else:              # duplication of single base
		try：
			locus = re.search(r'(\d+)dup', posedit)
			num_left =  int(locus.group(1)) 
			num_right = num_left
			num_bases = False 
		except:
			sys.exit(f'Error: parse {posedit} error')

	ref_from_string = posedit.split('dup')[1]
	ref_from_locus = seq[num_left - 1:num_right]
	if ref_from_string != ref_from_locus:
		sys.exit(f'Warning: {posedit} is different from seq[start - 1:end]', file=sys.stderr, flush=True)
	
	# before duplication
	seq_left10 = seq[num_left - 11:num_left - 1]
	seq_right10 = seq[num_right :num_left + 9]
	duplication = seq[num_left - 1:num_right]
	seq20_before_dup = seq_left10 + duplication + seq_right10
	# after duplication
	seq_after_dup = seq_left10 + duplication + duplication + seq_right10

	step = len(duplication)
	seq_len = len(seq20__before_dup)	
	for i in range(10 + step, seq_len, step):
		shift_dup = seq[i: i + step]
		if shift_dup != duplication:
			break
		shift = True

	locus = [num_left + i - 10, num_left + i - 10 + step - 1]
	if shift:
		if num_bases:
			print(f'old: {num_left}_{num_right} -> new: {locus[0]}_{locus[1]}', file=sys.stdout, flush=True)
		else:
			print(f'old: {num_left} -> new: {locus[0]}', file=sys.stdout, flush=True)
	
	if num_bases:
		return ''.join(['_'.join([str(pos) for pos in locus]), 'dup', ref_from_locus])
	else:
		return  ''.join([str(locus[0]), 'dup', ref_from_locus])

def priority_dup_lt_ins(seq:str, ins_posedit:str) -> str:
	# deletion > duplication > insertion > deletion–insertion.
	if '_' in posedit:
		locus = re.seqrch(r'(\d+)_(\d+)ins', ins_posedit)
		num_left = locus.group(1)
		num_right = locus.group(2)
		num_bases = True
	else:
		locus = re.search(r'(\d+)ins*', ins_posedit)
		num_left = locus.group(1)
		num_right = num_left
		num_bases = False
	else:
		pass

	ref_from_string = posedit.split('ins')[1].strip('\n').strip()
	length_ref_from_string = len(ref_from_string) 
	
	ref_left_bases = seq[num_left - length_ref_from_string:num_left]
	if ref_left_bases == ref_from_string:
		change = True
		locus = [num_left - length_ref_from_string + 1, num_left]
		if num_bases:
			return ''.join(['_'.join([str(pos) for pos in locus]), 'dup', ref_from_string])
		else:
			return ''.join([str(locus[0]), 'dup', ref_from_string])
	else:
		return posedit
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys
import  difflib
import argparse
import os
import pandas as pd


def correct_orientations(sequences_fasta, filename, forward, reverse):
	forward_comp = str(Seq(forward).reverse_complement())
	reverse_comp = str(Seq(reverse).reverse_complement())


	sequences2save = []
	for pos, record in enumerate(SeqIO.parse(sequences_fasta, "fasta")):
		seq = record.seq
		
		# forward -> reverse
		left_primer_fr = seq[:len(forward)]
		right_primer_fr = seq[-len(reverse):]

		
		matches_left_fr = difflib.get_close_matches(left_primer_fr, [forward, forward_comp, reverse, reverse_comp])
		best_match_left_fr = matches_left_fr[0]

		matches_right_fr = difflib.get_close_matches(right_primer_fr, [forward, forward_comp, reverse, reverse_comp])
		best_match_right_fr = matches_right_fr[0]

		# reverse -> forward
		left_primer_rf = seq[:len(reverse)]
		right_primer_rf  = seq[-len(forward):]

		
		matches_left_rf = difflib.get_close_matches(left_primer_rf, [forward, forward_comp, reverse, reverse_comp])
		best_match_left_rf = matches_left_rf[0]

		matches_right_rf = difflib.get_close_matches(right_primer_rf, [forward, forward_comp, reverse, reverse_comp])
		best_match_right_rf = matches_right_rf[0]


		# comparison of different orientations matches
		left_fr_score = difflib.SequenceMatcher(None, left_primer_fr, best_match_left_fr).ratio()
		right_fr_score = difflib.SequenceMatcher(None, right_primer_fr, best_match_right_fr).ratio()
		left_rf_score = difflib.SequenceMatcher(None, left_primer_rf, best_match_left_rf).ratio()
		right_rf_score = difflib.SequenceMatcher(None, right_primer_rf, best_match_right_rf).ratio()


		left_match = best_match_left_fr if left_fr_score > left_rf_score else best_match_left_rf
		right_match = best_match_right_fr if right_fr_score > right_rf_score else best_match_right_rf
		
		# check orientations

		record = SeqRecord(
			record.seq,
			id=str(pos + 1),
			description=''
		)
		f_r = False
		f_rc = False
		fc_r = False
		r_fc = False
		if left_match == forward and right_match == reverse:
			f_r = True
			record = SeqRecord(
	    		record.seq,
	    		id=str(pos + 1),
	    		description=''
			)
			sequences2save.append(record)


		elif left_match == forward and right_match == reverse_comp:
			f_rc = True
			reverse_complement = Seq(str(record.seq)).reverse_complement()
			record = SeqRecord(
	    		reverse_complement,
	    		id=str(pos + 1) + 'f_rc',
	    		description=''
			)
			sequences2save.append(record)

		elif left_match == forward_comp and right_match == reverse:
			fc_r = True
			reverse_complement = Seq(str(record.seq)).reverse_complement()

			record = SeqRecord(
	    		reverse_complement,
	    		id=str(pos + 1) + 'fc_r',
	    		description=''
			)
			sequences2save.append(record)

		elif left_match == reverse and right_match == forward_comp:
			r_fc = True
			record = SeqRecord(
	    		record.seq,
	    		id=str(pos + 1) + 'r_fc',
	    		description=''
			)
			sequences2save.append(record)

	print(filename, f_r, f_rc, fc_r, r_fc)

	SeqIO.write(sequences2save, f'corrected_orientations/{filename}', "fasta")


primers = pd.read_csv('primers_files/primers.tsv', sep='\t')
for file in os.listdir('fasta_files'):
	primer_name = file.replace('.fa', '')
	primer = primers.loc[primers['NAME'] == primer_name].reset_index(drop=True)
	forward = primer['FORWARD'][0]
	reverse = primer['REVERSE'][0]

	correct_orientations(f'fasta_files/{file}', file, forward, reverse)

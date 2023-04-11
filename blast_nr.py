'''

Written by Elinor 6/30/22 to grab all sequences in R2G files, write to fasta, BLASTx against NR, then write output to a tsv
input is gc3_master_length_and_complexity.csv created by get_composition_and_complexity.py
First part is new, then combines parseXML.py and BLASTseqsNR.py at the end.
this version takes into account a length and coverage cutoff as well

Updated 7/27/22

Goal is to specify and minimize which sequences to BLAST, since this takes a long time and is CPU intensive

Usage: Coverage and length are optional Run python3 grab_outliers.py for help

'''

import os
import sys
from Bio.Blast.NCBIWWW import qblast
from Bio import SeqIO
from tqdm import tqdm
from Bio.Blast import NCBIXML
from pathlib import Path



def script_help():
	print(f'This was written to extract specific sequences from ready to go files check how many outliers there are, and BLAST them. You tell it the ranges and cutoffs that you want (inliers) and it give you the outliers\n\n\t You can specify: \n\n --len <number> to only get sequences above a certain length. The default is 100bp \n\n --cov <number> to only get sequences above a certain coverage. The default is 1\n\n\t You must specify:\n\n a range of GC3 with --gc3_range <lower cutoff> <higher cutoff> like this: \n\n --gc3_range 3 13.3 \n\n and your clade of interest with --clade <1 to 10 digit code> \nlike this: \n\n --clade Sr_r\n')
	exit()

fasta_handle = 'ReadyToGo_NTD'


def get_args():

	blast = input('Do you want to BLAST your outliers? (y or n): ')
	if blast == 'y':
		BLASTseqsNR()
	else:
		exit()


def BLASTseqsNR():

	print('\n\nNow BLASTing sequences in to_BLAST_by_taxon. It will BLAST all sequences one taxon at a time. \n\nThis may take a while ...   BLAST may penalize you if you send too much :((')

	#Blasts each fasta file against the NR database and returns and one XML file of hits per query.
	Path(f'BLAST_records').mkdir(parents=True, exist_ok=True)


	#BLAST outliers by taxon (files in 'to_BLAST_by_taxon'), which sends fewer seqs per job. they do not run in parallel
	for file in os.listdir('LSGF_NTD'):
		if file.endswith('fa'):
			all_seqs = list(SeqIO.parse(f'LSGF_NTD/{file}', 'fasta'))
			print(f'\nBLAST-ing Sequences in {file} \n')
			for record in tqdm([s for s in all_seqs]):
				with open(f'BLAST_records/{record.id}', 'w') as o:
					results = qblast(program = 'blastx', sequence = record.seq, database = 'nr', hitlist_size = 8, format_type = 'XML')
					o.write(results.read())
					results.close


	parseXML()



def parseXML():

	print('\n\n\nNow parsing BLAST output to a TSV file \n\n\n')
	

	with open('BLASTnrOutliers.tsv', 'w') as o:
		o.write("Seq Title\tQuery Name\tQuery Length\tAlignment Title\tAlignment ID\tAlignment Def\teValue\tCoverage\tPercent_ID")
		for xml_file in os.listdir('BLAST_records'):
			if xml_file.startswith(clade):
				result_handle = open('BLAST_records/' + xml_file)
				blast_records = NCBIXML.parse(result_handle)
				for rec in blast_records:
					for alignment in rec.alignments:
						for hsp in alignment.hsps:
							o.write('\n'+str(xml_file) +'\t'+ str(rec.query_id) + '\t' + str(rec.query_length) + '\t' + str(alignment.title) + '\t' + str(alignment.hit_id) + '\t' + str(alignment.hit_def) + '\t' + str(hsp.expect) + '\t' + str(hsp.align_length / rec.query_length) + '\t' + str(hsp.identities/ hsp.align_length))

if __name__ == '__main__':
	if len(sys.argv) < 1:
		script_help()
	else:
		get_args()

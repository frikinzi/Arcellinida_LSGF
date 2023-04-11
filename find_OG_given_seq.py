import os
from Bio import SeqIO

for LSGF in os.listdir("LSGFs"):
    seqs = list(SeqIO.parse("LSGFs/" + LSGF, "fasta"))
    for seq in seqs:
        if os.path.exists("BLAST_records/" + seq.id):
            os.system("mv BLAST_records/" + seq.id + " BLAST_records_renamed/" + seq.id + "_" + LSGF)



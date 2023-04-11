import os
from Bio import SeqIO

newfasta = open("BLASTed.fa", "w")
for LSGF in os.listdir("LSGF_NTD"):
    if LSGF != ".DS_Store":
        seqs = list(SeqIO.parse("LSGF_NTD/" + LSGF, "fasta"))
        for seq in seqs:
            if os.path.exists("BLAST_records_renamed/" + seq.id + "_" + LSGF):
                newfasta.write(">" + LSGF + "_" + seq.description + "\n")
                newfasta.write(str(seq.seq) + "\n\n")



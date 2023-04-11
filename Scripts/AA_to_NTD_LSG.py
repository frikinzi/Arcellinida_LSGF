"""Given two folders of AA and NTD sequences, matches the AA sequence to the NTD sequence, and writes the NTD sequences to a new fasta file

Author: Angela Jiang"""

import os
from Bio import SeqIO

big_dict = {}
for taxon in os.listdir("All_TIdeS_ORFs"):
    if taxon.endswith(".fas"):
        seqs = list(SeqIO.parse("All_TIdeS_ORFs/" + taxon, "fasta"))
        for seq in seqs:
            big_dict["_".join(seq.description.split("_")[:7])] = seq.seq


for LSGF in os.listdir("LSGF_Sequences/LSGFs_elegans/"):
    if LSGF != ".DS_Store":
        newfasta = open("LSGFs_elegans_NTD/"+ LSGF, "w")
        seqs = list(SeqIO.parse("LSGF_Sequences/LSGFs_elegans/" + LSGF, "fasta"))
        for seq in seqs:
            try:
                newfasta.write(">" + seq.description + "\n")
                newfasta.write(str(big_dict["_".join(seq.description.split("_")[:7])]) + "\n\n")
            except:
                pass


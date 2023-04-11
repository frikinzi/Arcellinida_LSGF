"""Given a folder of fasta sequences, creates a gene count table like the output of OrthoFinder. Put a list of your taxa in a text file called taxon_list.txt

Usage: python3 creategenecounttable.py [NAME OF FOLDER] """

import os
from Bio import SeqIO
import sys

taxon_list = open("taxon_list.txt", "r")

listtaxa = taxon_list.readlines()

taxonlist = []

for a in listtaxa:
    if len(a) > 8:
        taxonlist.append(a[:10].strip())

taxon_list.close()

newtable = open("curated_genecount.csv", "w")
towrite = []
for taxa in taxonlist:
    towrite.append(taxa + ",")
newtable.write("Orthogroup,")
for item in towrite:
    newtable.write(item)
newtable.write("\n")

listofogs = os.listdir(sys.argv[1])

for og in listofogs:
    if og.endswith("fasta") or og.endswith("fa"):
        newtable.write(og[:-20] + ",")
        #create dictionary with taxon:0
        count_dict = {}
        for taxon in taxonlist:
            count_dict[taxon] = 0
        seqlist = list(SeqIO.parse(sys.argv[1]+ "/" + og, "fasta"))
        for seq in seqlist:
            try:
                count_dict[seq.description[:10]] = count_dict.get(seq.description[:10]) + 1
            except:
                pass
        for key,value in count_dict.items():
            newtable.write(str(value) + ",")
        newtable.write("\n")

        
            


    

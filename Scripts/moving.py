"""Given a folder of fasta sequences and a list of OGs in a text file called OG_list.txt, moves them to a new folder

Usage: python3 moving.py [NAME OF FOLDER] [NAME OF FOLDER YOU WANT THEM MOVED TO]"""

import os
import sys

og_file = open("OG_list.txt", "r")
og_list = og_file.readlines()

for i in og_list:
    os.system("cp " + sys.argv[1] + "/" + i.strip() + "aln.With_Names.fasta " + sys.argv[2] + "/")
"""Given a folder of sequences (folder named curatedround0), removes sequences more than 2x the median, and 0.5x the median, and iterates the curation 11 times

Author: Angela Jiang"""

import os
from Bio import SeqIO
import statistics


for y in range(11):
    listdir = os.listdir("curatedround" + str(y))
    try:
        os.mkdir("curatedround" + str(int(y)+1))
    except:
        pass

    for i in listdir:
        if i.endswith("fasta") or i.endswith("fa"):
            newfile = open("curatedround" + str(int(y)+1) + "/" + i, "w")
            seqlist = list(SeqIO.parse("curatedround" + str(y) + "/" + i, "fasta"))
            listoflengths = []
            for a in seqlist:
                listoflengths.append(len(str(a.seq)))
            try:
                average = statistics.median(listoflengths)
                for a in seqlist:
                    if not (len(str(a.seq)) < 0.5 * average or len(str(a.seq)) > 2 * average):
                        newfile.write(">" + a.description + "\n")
                        newfile.write(str(a.seq) + "\n" + "\n")
            except:
                print("there are 0 sequences in " + i)
                pass
            newfile.close()
        

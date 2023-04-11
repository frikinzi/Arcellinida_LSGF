"""A wrapper of iqtree. runs iqtree on a folder of sequences

Author: Angela Jiang"""

import os
import sys

for i in os.listdir(sys.argv[1]):
    if (i.startswith("OG")):
        os.system("iqtree2 -s " + sys.argv[1]+ "/" + i + "  -B 1000 -bnni -ntmax 12 -nt AUTO")
'''
Plots a panel of codon usage plots given tables of codon usage data

Author: Angela Jiang
'''

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def main():
	listofdir = os.listdir("codon_usage_folder")
	listofdir.sort()
	inte = 1
	fig = plt.figure()
	#plt.subplots(nrows=2, ncols=5)
	for a in listofdir:
		if (a != ".DS_Store"):
			plt.subplot(5, 8, inte)
			null = pd.read_table("codon_usage_folder/" + a + "/SpreadSheets/" + a[:12] + ".ENc.Null.tsv")
			pz = pd.read_table("codon_usage_folder/" + a + "/SpreadSheets/" + a[:12] + ".ENc.Raw.tsv")
			df1, df2 = pz.loc[pz['LSG'] == 1 ] , pz.loc[pz['LSG'] == 0 ]
			sns.lineplot(data=null, x='GC3', y='ENc', lw=2,color='black')
			sns.scatterplot(data=df1, x='GC3-Degen', y='ObsWrightENc_6Fold', s=12, color = "orange")
			sns.scatterplot(data=df2, x='GC3-Degen', y='ObsWrightENc_6Fold', s=12, color = "blue")
			plt.xlim(0,101)
			plt.ylim(20,66)
			
			plt.xlabel("")
			plt.ylabel("")
			plt.xticks(fontsize=8)
			plt.yticks(fontsize=8)
			plt.title(a[:12], fontsize=8)
			inte = inte + 1
		
	plt.tight_layout(pad=0.3, w_pad=-1, h_pad=-1)
	fig.supxlabel('GC3')
	fig.supylabel('ENc')
	fig.legend()
	plt.show()
	fig.savefig("jellyfishplotpanels.png")

main()

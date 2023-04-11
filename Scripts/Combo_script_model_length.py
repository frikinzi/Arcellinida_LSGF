"""Intended to help with basic plotting (box/violinplots) in python.

Seaborn boxplot: https://seaborn.pydata.org/generated/seaborn.boxplot.html

Seaborn violinplot: https://seaborn.pydata.org/generated/seaborn.violinplot.html

Above the links to help. They have great examples if needed!"""
#Script by Emma Schumacher
#Last updated/compiled by ETS 05/16/22

#Import statements...
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from random import sample
import sys
from pathlib import Path
import os
import shutil

# my main function
def main(og_fasta_file, og_number, random, run, multi, s_size):
    # maybe condense randoms ???
    if (random == True):
        # if only one run and not more than one fasta, does randomly sample and goesnt call make_many_OG_boxplots
        if ((run == 1) & (multi == False)): 
            functionscall(og_fasta_file, og_number, boxplot_with_seaborn, multi, random, samp = s_size)            
        else:
            functionscall(og_fasta_file, og_number, make_many_OG_boxplots, multi, random, samp = s_size, runs = run)
    # if not random 
    else:
        print("not")
        if ((run == 1) & (multi == False)):
            functionscall(og_fasta_file, og_number, boxplot_with_seaborn, multi, random)
        else:
            functionscall(og_fasta_file, og_number, make_many_OG_boxplots, multi, random, runs = run)

def functionscall(file, og_number, graphprep, multi, random, samp = None, runs = None):
    #if there is random sampling
    if (random):
        #if more than one fasta
        if (multi):
            #Randomly samples from all files
            some_fasta = []
            for fasta in file:
                #print("ran init")
                print(fasta)
                a_fasta = make_multiple_randoms(fasta, samp, runs)
                some_fasta.append(a_fasta)
                
        #if only one fasta
        else:
            #Randomly samples from the single file 
            some_fasta = make_multiple_randoms(file, samp, runs)
            
    #if there is no scrambling, doesn't do anything to some_fasta
    else:
        some_fasta = file
        
    
    #Calls function to make a tsv for next function
    #print("ran tsv")
    some_tsv = make_multiple_OGs_sizes(some_fasta, og_number, runs, random)
    #Makes a boxplot
    #print("ran graph")
    print(some_tsv)
    graphprep(some_tsv) 

#used to randomly sample 500 OGs from fasta to make stuff run better
def get_random(og_fasta_file, samp, diff = None):
    path = os.getcwd() 

    #if there is only one run
    if not diff:
        print(og_fasta_file)
        corrected_file = path + og_fasta_file[og_fasta_file.find('/', len(path) + 1):og_fasta_file.index(".fasta")] + "smallerseq.fasta"
        #corrected_file = og_fasta_file[0:og_fasta_file.index(".fasta")] + "smallerseq.fasta"
        
    #if multiple runs, changes names of files based on iteration
    else:
        corrected_file = path + og_fasta_file[og_fasta_file.find('/', len(path) + 1):og_fasta_file.index(".fasta")] + str(diff) + "smallerseq.fasta" 

    
    #opens original file and new file to write 500 samples from the latter to the former
    with open(og_fasta_file, "r") as original, open(corrected_file, 'w') as corrected:
        
        seqs = SeqIO.parse(og_fasta_file, "fasta")
        
        for seq in sample(list(seqs), samp): #no replacement ??? too big???
            SeqIO.write(seq, corrected, 'fasta')
            
    return(corrected_file)
    # ??? PROPORTION PRESERVATION

# Provide a FASTA file of sequences, along with an OG number, to make a quick
# table of sequence lengths (note this is per OG).
def get_sizes(fasta_file, og_number):
    
    # Snag the sequence names and their lengths from a given FASTA file. 
    seq_sizes = {i.id:len(i.seq) for i in SeqIO.parse(fasta_file,'fasta')}
    
    # Save the data as a spreadsheet with Tab-Separated-Values (TSV)
    with open(f'{og_number}.SeqLength.tsv','w+') as w:
        w.write('OG\tSequence_Name\tLength\tDomain\n')

        for k, v in seq_sizes.items():
            # Add in the taxonomic "domains" when possible! This is mostly for
            # fancy plotting. 
            if k[:2] == 'Ba':
                w.write(f'{og_number}\t{k}\t{v}\tBacteria\n')

            else:
                w.write(f'{og_number}\t{k}\t{v}\tEukaryota\n')
        return(og_number + ".SeqLength.tsv")

#used to call size tsvs for multiple files 
def make_multiple_OGs_sizes(fasta_files, og_number, runs, random = True):
    #list of all the tsvs we make
    all_OG_tsvs = []

    #if random and not a curated list, changes naming
    if (random): 
        if (isinstance(og_number, list)):
            for fasta in fasta_files:
                og = og_number.pop()
                if (runs > 1):
                    for i in range(0, runs): 
                        fn = get_sizes(fasta.pop(), (og + "_sample_" + str(i))) 
                        all_OG_tsvs.append(fn)
                else: 
                    fn = get_sizes(fasta, og)
                    all_OG_tsvs.append(fn)
                    
        else:
            if ((runs == 1) | (runs == None)):
                fn = get_sizes(fasta_files, og_number) 
                all_OG_tsvs.append(fn)
            else:
                for i in range(0, runs):
                    fn = get_sizes(fasta_files.pop(), (og_number + "_sample_" + str(i)))
                    all_OG_tsvs.append(fn)

    #else, names based on curated list you pass in
    else:
        if (isinstance(og_number, list)):
            for fasta in fasta_files:
                og = og_number.pop()
                fn = get_sizes(fasta, og) #pop?
                all_OG_tsvs.append(fn)
        else: 
            fn = get_sizes(fasta_files, og_number)
            all_OG_tsvs.append(fn)
    print(all_OG_tsvs)
    return (all_OG_tsvs)

#used to call multiple random samples, great for seeing if distributions are consistent
def make_multiple_randoms(file, samp, runs):
    if ((runs == 1) | (runs == None)):
        all_samples = get_random(file)
    else:
        #holds all the random samples
        all_samples = []

        #gets that many randoms (with resampling)
        for i in range(0, runs):
            fn = get_random(file, samp, diff = i)
            all_samples.append(fn)
    return (all_samples)

    #??? GET OUTLIERS

# Generate a boxplot with seaborn, for many OGs at once 
def make_many_OG_boxplots(all_OG_tsvs):
    #reads first tsv, doesnt remove header line 
    f = pd.read_csv(all_OG_tsvs.pop(0), sep = '\t', header = 0)
    #creates list of tsv pandas dataframes
    newfiles = [f]

    #for all the tsvs we made other than the first one
    for tsv in all_OG_tsvs: 
        #turn tsv into dataframe
        f = pd.read_csv(tsv, sep = '\t', header = 0)
        #removes header so it isn't repeated a million times
        f = f.iloc[:-1]
        #adds new dataframe to list of pandas dataframes we want to combine
        newfiles.append(f)
    
    #combines all the tsv dataframes
    f = pd.concat(newfiles)
    
    #makes an actual file
    f.to_csv('combined.tsv', sep="\t")

    #calls boxplot
    boxplot_with_seaborn('combined.tsv') 
    

# A couple examples of different approaches/options for generating boxplots and
# violinplots in python.
def boxplot_with_seaborn(some_tsv): 
    slen_df = pd.read_table(some_tsv, engine = 'python')
    #slen_df = pd.read_csv(some_tsv, sep = '\t', header = 0, engine = 'python') #???
    #gets graph axis
    
    #uprlim = slen_df['Length'].max() + 50 
    column = slen_df["Length"] 
    uprlim = column.max() + 50
    print(uprlim)
    #lowlim = slen_df["Length"].min()

    # Plot the boxplot to visualize the lengths (essentially "raw")
    ax = sns.boxplot(data=slen_df, x='OG', y='Length')

    # Adjust the y-axis (min, max)
    ax.set_ylim(0, uprlim)

    # Show the plot. You can then save it if you want, otherwise just close it!
    plt.show()

    # Plot the boxplot to visualize the lengths by taxonomic domain!
    ax = sns.boxplot(data=slen_df, x='OG', y='Length', hue='Domain')
    plt.show()

    # Plot the boxplot to visualize the lengths by taxonomic domain, without outliers.
    ax = sns.boxplot(data=slen_df, x='OG', y='Length', hue='Domain', showfliers=False)
    plt.show()

    # Violinplot to visualize the lengths by taxonomic domain, without outliers.
    ax = sns.violinplot(data=slen_df, x='OG', y='Length', showfliers=False)
    #trying to add ticks
    #plt.yscale('log')
    #plt.yticks(range(0, 2500, 100))
    plt.show()

def deal_w_input(path, og_number, s_size = 0, runs = 0):
    # gets input 
    og_fasta_file = []
    for file in os.listdir(path):
        if not file.startswith('.'):
            file = path + file
            og_fasta_file.append(file)
    
    # checks if multiple folder
    if (len(og_number) > 1):
        multi = True
    else:
        multi = False

    #checks for obvious mistake
    if (len(og_fasta_file) != len(og_number)):
        print("You need an og number for each fasta file/you need a fasta file for each og number")
    if (len(og_number) == 1):
        og_fasta_file = str(og_fasta_file.pop())
        og_number = str(og_number.pop())
    
    #if you want samples ???
    if ((s_size > 0) & (runs > 0)):
        random = True
    else:
        random = False
    
    #calls my main function
    main(og_fasta_file, og_number, random, runs, multi, s_size)

 #if len(l)<4 else random.sample(l, 4)??? 
if __name__ == '__main__':
    ## EDIT HERE: Put your fasta files into an input folder and add the path
    path = sys.argv[1] + "/"
    
    ## EDIT HERE: Modify this so it is whatever label you want on your X axis, if you want to run
    # multiple different use a list of names []
    listofdir = os.listdir(path)
    og_number = []
    for a in listofdir:
    	if a.endswith("fa"):
        	og_number.append(a[:-3])
    
    

    ## EDIT HERE: Modify this so it is whatever sample size you want if you want random samples,
    # if you dont want random samples, make this 0
    s_size = 0
    
    ## EDIT HERE: Modify this so it is however many runs you want (only matters if
    # you want random batches
    runs = 0
    
    # deals with input
    deal_w_input(path, og_number, s_size, runs) 






    

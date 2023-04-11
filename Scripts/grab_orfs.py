#!/usr/bin/env python3

''' This script grabs the largest ORF from a given transcript, using one of
several valid genetic codes. User-defined minimum lengths also cull the putative
set of ORFs.

This is similar to "TransDecoder", which did not support alternative genetic
codes, which are prevalent across the eukaryotic tree of life.'''

# Author: Xyrus Maurer-Alcala
# Contact: maurerax@gmail.com
# Last Modified: 2020-07-20
# usage: python

# Dependencies:
# Python3, Barrnap, CD-HIT, BioPython

import logging
import os
import re
import shutil
import subprocess
import sys


from collections import defaultdict
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def check_dependencies():
    d = {"Barrnap": shutil.which("barrnap"),
        "CD-HIT-EST":shutil.which("cd-hit-est")}
    return d


def prep_folders(taxon_code):
    taxon_dir = f'{taxon_code}_ORF_Call'

    Path(f'{taxon_dir}/Clustered_Transcripts').mkdir(parents=True, exist_ok=True)
    Path(f'{taxon_dir}/rDNA_Seqs').mkdir(parents=True, exist_ok=True)
    Path(f'{taxon_dir}/ORFs').mkdir(parents=True, exist_ok=True)
    Path(f'{taxon_dir}/ORFs').mkdir(parents=True, exist_ok=True)
    Path("All_ORFs/NTD").mkdir(parents=True, exist_ok=True)
    Path("All_ORFs/AA").mkdir(parents=True, exist_ok=True)

    return taxon_dir


def check_trans_table(gcode):
    valid_trans_table = {1:["TGA","TAG","TAA"], 4:["TAG","TAA"], 6:["TGA"],
        10:["TAA","TAG"], 12:["TGA","TAG","TAA"], 29:["TGA"], 30:["TGA"],
        "chilo":["TAA"]}

    if gcode in valid_trans_table.keys():
        return valid_trans_table[gcode]

    else:
        print("Invalid genetic code provided. Defaulting to the universal " \
            "tranlsation table.")
        return 1


# If a log/run has been stared already on the same day, then update the
# suffix to ensure subsequent runs do not overwrite current logs files.
def check_log_folder_exist(path):
    filename, extension = os.path.splitext(path)
    counter = 2

    while os.path.exists(path):
        path = f'{filename}.Run_{counter}{extension}'
        counter += 1

    return path


def param_depend_log(input_trans, taxon_code, gcode, min_orf):

    curr_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")

    taxon_dir = prep_folders(taxon_code)

    # Initializes the log with the date and increments the number of the suffix
    # if a run has already been initiated on the same day.
    current_date_log = check_log_folder_exist(
        f'{taxon_dir}/{taxon_code}.{datetime.now().strftime("%b_%d")}.log')

    logging.basicConfig(filename=current_date_log,
        format='%(message)s',
        level=logging.DEBUG)

    # Saves the relevant information to the log file.
    logging.info(f'Thanks for trying to grab some ORFs!')
    logging.info(f'\nORF-calling run, started: {curr_time}')

    logging.info(
        f'\n+---------------------------------+\n'\
        f'|       ORF-Call Parameters       |\n'\
        f'+---------------------------------+'
        )

    trans_name = input_trans.split("/")[-1]

    logging.info(f'   Transcriptome: {trans_name}')
    logging.info(f'   Translation table: {gcode}')
    logging.info(f'   Taxon Code: {taxon_code}')
    logging.info(f'   ORF Size (bp): {min_orf}')

    dependencies = check_dependencies()

    logging.info(
        f'\n+---------------------------------+\n'\
        f'|      Checking Dependencies      |\n'\
        f'+---------------------------------+'
        )

    for k, v in dependencies.items():
        m = 0
        if v:
            logging.info(f'   {k}: check')
        else:
            logging.warning(f'   {k}: missing')
            m += 1
        if m > 0:
            logging.warning(f'    Aborting run given missing dependencies.')
            sys.exit(1)


    logging.info(
        f'\n+---------------------------------+\n'\
        f'|   Checking Translation Table    |\n'\
        f'+---------------------------------+'
        )

    if gcode.isdigit():
        trans_table = check_trans_table(int(gcode))
    else:
        gcode = gcode.lower()
        trans_table = check_trans_table(gcode)

    logging.info(f'   Stop Codons: {", ".join(trans_table)}')

    return taxon_dir, trans_table


def call_rDNA(clustered_trans, taxon_code, out_dir):

    print(
        f'\n+---------------------------------+\n'\
        f'|     Removing Putative rRNAs     |\n'\
        f'+---------------------------------+'
        )

    logging.info(
        f'\n+---------------------------------+\n'\
        f'|     Removing Putative rRNAs     |\n'\
        f'+---------------------------------+'
        )

    euk_out = f'{out_dir}/rDNA_Seqs/Euk_rRNA.fas'
    mito_out = f'{out_dir}/rDNA_Seqs/Mito_rRNA.fas'
    bac_out = f'{out_dir}/rDNA_Seqs/Bact_rRNA.fas'
    arc_out = f'{out_dir}/rDNA_Seqs/Arch_rRNA.fas'
    clean_out = f'{out_dir}/rDNA_Seqs/{taxon_code}.No_rDNA.fas'

    rDNA_seqs = []

    rDNA_db = {"bac":bac_out,"arc":arc_out, "mito":mito_out, "euk":euk_out}
    for k, v in rDNA_db.items():
        barrnap_cmd = f'barrnap --kingdom {k} --outseq {v} {clustered_trans}'

        barrnap_call = subprocess.call(barrnap_cmd, shell=True,
            stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    for v in rDNA_db.values():
        rDNA_seqs += [i.id.split(":")[2] for i in SeqIO.parse(v,"fasta")]
    logging.info(f'   Number of rRNA seqs: {len(set(rDNA_seqs))}')

    with open(clean_out, "w+") as w:
        for i in SeqIO.parse(clustered_trans,"fasta"):
            if i.description not in rDNA_seqs:
                w.write(f'>{i.description}\n{i.seq}\n')

    return clean_out


def cluster_transcripts(transcriptome, taxon_code, out_dir, round=1, minLen=300):
    if round == 1:
        shutil.copy(transcriptome, out_dir)

        print(
            f'\n+---------------------------------+\n'\
            f'|   "Raw" Transcript Clustering   |\n'\
            f'+---------------------------------+'
            )

        logging.info(
            f'\n+---------------------------------+\n'\
            f'|   "Raw" Transcript Clustering   |\n'\
            f'+---------------------------------+'
            )

        pre_clust_fas = f'{out_dir}/Clustered_Transcripts/{taxon_code}.ToClust.fas'
        post_clust_fas =f'{out_dir}/Clustered_Transcripts/{taxon_code}.CdHit.fas'
        pre_clust_codes = f'{out_dir}/Clustered_Transcripts/{taxon_code}' \
            f'.ToClust.SeqCodes.tsv'

        trans_seqs = [i for i in SeqIO.parse(transcriptome,"fasta") if len(i.seq) >= minLen]
        trans_seqs.sort(key=lambda x: -len(x.seq))

        seq_codes = {}
        seq_count = 1

        with open(pre_clust_fas,"w+") as w:
            for i in trans_seqs:
                w.write(f'>{taxon_code}_Trans_{seq_count}\n{i.seq}\n')
                seq_codes[i.id] = f'{taxon_code}_Trans_{seq_count}'
                seq_count += 1

        with open(pre_clust_codes, "w+") as w:
            w.write("Original_SeqName\tRenamed_Seqn\n")
            for k, v in seq_codes.items():
                w.write(f'{k}\t{v}\n')

        logging.info(f'    Pre-clustered transcripts: {len(trans_seqs)}')
        logging.info(f'    "Clean" transcripts: {seq_count-1}')

    else:
        pre_clust_fas = transcriptome
        post_clust_fas = transcriptome.replace("Initial","Final")

        print(
            f'\n+---------------------------------+\n'\
            f'|          ORF Clustering         |\n'\
            f'+---------------------------------+'
            )

        logging.info(
            f'\n+---------------------------------+\n'\
            f'|          ORF Clustering         |\n'\
            f'+---------------------------------+'
            )

    cd_hit_est_cmd = f'cd-hit-est -G 0 -c 0.97 -aS 1.0 -aL 0.005 -i ' \
        f'{pre_clust_fas} -o {post_clust_fas}'

    cd_hit_call = subprocess.call(cd_hit_est_cmd, shell=True,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    if round == 2:
        fin_orfs = len([i for i in SeqIO.parse(post_clust_fas,'fasta')])
        logging.info(f'    There are {fin_orfs} "clean" ORFs')

    return post_clust_fas


def extractORFs(stopCDNs, seq, minLen=300):
    '''Will return the single largest complete (start and stop codon inclusive)
    ORF present in a transcript that meets the minimum length criteria'''

    regex = re.compile(r'(' + "|".join(stopCDNs) + r')')
    inFrame_stops = list(set([0]+[m.end() for m in regex.finditer(seq) if m.end()%3 == 0]))
    inFrame_stops.sort()
    dirtyORFs = [seq[i:j] for i,j in zip(inFrame_stops, inFrame_stops[1:] + [None])]
    # Keeps the putative ORFs that meet and/or exceed the minimal protein length.
    # clean_ORFs = [i for i in ORFs if len(i) >= minLen*3]
    cleanORFs = []
    for i in dirtyORFs:
        if len(i) >= minLen and i[-3:] in stopCDNs:
            # start_pos = [m.start() for m in re.finditer("ATG", i[:60]) if m.start()%3 == 0]
            start_pos = [m.start() for m in re.finditer("ATG", i) if m.start()%3 == 0]
            if start_pos:
                withStart = i[start_pos[0]:]
                if len(withStart) >= minLen and len(withStart)%3 == 0:
                    cleanORFs += [withStart]

    if cleanORFs:
        return max(cleanORFs,key=len)
    else:
        return None


def collect_orfs(slim_fasta, stop_cdns, taxon_code, out_dir, min_orf = 300):

    print(
        f'\n+---------------------------------+\n'\
        f'|         Identifying ORFs        |\n'\
        f'+---------------------------------+'
        )

    logging.info(
        f'\n+---------------------------------+\n'\
        f'|         Identifying ORFs        |\n'\
        f'+---------------------------------+'
        )

    ORF_fasta = f'{out_dir}/ORFs/{taxon_code}.InitialORFs.fas'
    largest_ORFs = defaultdict(int)
    final_ORFs = {}

    for i in SeqIO.parse(slim_fasta, "fasta"):
        revseq = str(i.seq.reverse_complement())
        fwdseq = str(i.seq)
        rfDict = {"_RF1":fwdseq, "_RF2":fwdseq[1:], "_RF3":fwdseq[2:],
            "_RF4":revseq, "_RF5":revseq[1:], "_RF6":revseq[2:]}
        for k, v in rfDict.items():
            # Links the orfs from a given reading frame to its original sequence
            temp_ORF = extractORFs(stop_cdns, v, min_orf)
            if temp_ORF:
                if len(temp_ORF) > largest_ORFs[i.description]:
                    largest_ORFs[i.description] = len(temp_ORF)
                    final_ORFs[i.description] = f'>{i.description}_Len_{len(temp_ORF)}{k}\n{temp_ORF}'

    with open(ORF_fasta,"w+") as w:
        w.write("\n".join(final_ORFs.values()))

    logging.info(f'    Extracted {len(final_ORFs)} putative ORFs larger than {min_orf}')

    return ORF_fasta


def back_up_ORFs(final_fasta, gcode):
    final_aa_fasta = final_fasta.replace(".fas",".AA.fas")

    with open(final_aa_fasta,"w+") as w:
        for i in SeqIO.parse(final_fasta,"fasta"):
            if gcode == "chilo":
                prot_seq = f'{i.seq[:-3].translate()}'.replace("*","Q")
            else:
                prot_seq = f'{i.seq[:-3].translate(gcode)}'
            w.write(f'>{i.description}\n{prot_seq}\n')

    shutil.copy(final_fasta, "All_ORFs/NTD/")
    shutil.copy(final_aa_fasta, "All_ORFs/AA/")

if __name__ == "__main__":
    if len(sys.argv)  < 4:
        print("\nUsage:")
        print("\ngrabORFs.py MyFasta MyTaxon genetic_code some_length")
        print("""\nNote: "some_length" refers to the DNA length, """\
            "so ensure it is a multiple of 3!")
        print("\nGenetic Codes:\n")
        gcd = ["1", "4", "6", "10", "12", "29", "30", "chilo"]
        print("\n".join(gcd)+"\n")
        sys.exit(0)

    args = sys.argv[1:]

    input_trans = args[0]
    taxon_code = args[1]
    gcode = args[2]

    if len(args) == 4:
        min_orf = max(6,int(args[3]))
        if min_orf%3 != 0:
            print("Minimum ORF size is not a multiple of 3. Aborting run!")
            sys.exit()
    else:
        min_orf = 300



    out_dir, trans_table = param_depend_log(input_trans, taxon_code, gcode, min_orf)


    initial_clust_fasta = cluster_transcripts(input_trans, taxon_code, out_dir)
    rDNA_cleaner = call_rDNA(initial_clust_fasta, taxon_code, out_dir)
    second_clust_fasta = collect_orfs(rDNA_cleaner, trans_table, taxon_code,
                            out_dir, min_orf)

    final_ORF_fasta = cluster_transcripts(second_clust_fasta, taxon_code, out_dir,2)

    back_up_ORFs(final_ORF_fasta, gcode)

    curr_time = datetime.now().strftime("%m/%d/%Y %I:%M %p")
    logging.info(f'\n     Finished! {curr_time}')
    print(f'Finished! {curr_time}')

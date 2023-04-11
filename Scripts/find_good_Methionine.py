"""Given a folder of alignments files, finds the first methionine that is shared by 50% of the sequences, then truncates any sequence that is more than 1.5x the median length to only include the methionine and the sequence after

Author: Angela Jiang"""

from Bio import AlignIO
import os
import csv


dir = (x for x in os.listdir("post_guidance_papilio_alignments/Alignments_post_guidance") if x.endswith("_Names"))

for file in dir:
    # Load the alignment from a FASTA file
    alignment = AlignIO.read("post_guidance_papilio_alignments/Alignments_post_guidance/" + file, "fasta")

    # Get the length of the alignment
    alignment_length = alignment.get_alignment_length()

    # Initialize a dictionary to store the counts of M at each position
    m_counts = {}

    # Loop over each position in the alignment
    for i in range(alignment_length):
        # Initialize a dictionary to store the counts of each amino acid at this position
        counts = {}
        # Loop over each sequence in the alignment
        for record in alignment:
            # Get the amino acid at this position in this sequence
            amino_acid = record.seq[i]
            # Increment the count for this amino acid
            if amino_acid in counts:
                counts[amino_acid] += 1
            else:
                counts[amino_acid] = 1
        # Get the count of M at this position
        if "M" in counts:
            m_counts[i] = counts["M"]

    #most_shared_position = max(m_counts, key=m_counts.get) just takes most shared pos

    min_m_count = 0.5 * len(alignment)
    most_shared_position = None
    temp_shared = max(m_counts, key=m_counts.get)
    for position in sorted(m_counts.keys()):
        if m_counts[position] >= min_m_count:
            print("found M position with more than 20 percent sequences " + str(position) + ", compared to max value " + str(temp_shared))
            most_shared_position = position
            break
    if most_shared_position == None:
        print("No more than 20 percent shared methionine position found in " + file)
        most_shared_position = temp_shared

    if temp_shared < most_shared_position:
        most_shared_position = temp_shared

    print(str(most_shared_position) + " for " + file)

    # Find the median length of the sequences
    lengths = []
    for record in alignment:
        length = sum(1 for aa in record.seq if aa != "-")
        lengths.append(length)
    median_length = sorted(lengths)[len(lengths)//2]

    # Initialize a list to store information about truncated sequences
    truncated_sequences = []

    # Loop over each sequence in the alignment
    new_records = []
    for record in alignment:
        # Check if the sequence is 1.5x the median length
        #print(str(len(record.seq))+ " median: " + str(median_length))
        if sum(1 for aa in record.seq if aa != "-") < 1.5 * median_length:
            new_records.append(record)
            continue
        # Find the nearest M to the most shared M
        nearest_m_position = None
        nearest_m_distance = float("inf")
        for i in range(most_shared_position - 30, most_shared_position + 31):
            if i < 0 or i >= alignment_length:
                continue
            if "M" == record[i]:
                distance = abs(i - most_shared_position)
                if distance < nearest_m_distance:
                    nearest_m_position = i
                    nearest_m_distance = distance
        # Truncate the sequence if an M is found, delete it otherwise
        if nearest_m_position is not None:
            new_seq = record.seq[nearest_m_position:]
            new_record = record[:]
            new_record.seq = new_seq
            if (sum(1 for aa in new_record.seq if aa != "-") >= 100):
                new_records.append(new_record)
                if (nearest_m_distance != 0):
                    truncated_sequences.append((record.id, "closest methionine at " + str(nearest_m_position) + " shared: " + str(most_shared_position)) )
            else:
                truncated_sequences.append((record.id, "too short after truncating"))
        else:
            truncated_sequences.append((record.id, "no_M_found"))
            continue

    # Create a new alignment with the remaining sequences
    new_alignment = alignment[:]
    new_alignment._records = new_records

    # Write the new alignment to a file
    AlignIO.write(new_alignment, "post_guidance_papilio_alignments/Curated/" + file + ".fasta", "fasta")

    # Write the summary CSV file
    with open("post_guidance_papilio_alignments/Tables_post/"+ file + ".csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Sequence ID", "Reason for truncation"])
        writer.writerows(truncated_sequences)

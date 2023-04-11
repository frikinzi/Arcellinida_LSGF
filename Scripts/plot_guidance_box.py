"""Given three concatenated Guidance score tables, creates a box plot of Guidance scores

Author: Angela Jiang"""

import pandas as pd
import matplotlib.pyplot as plt

# Set file paths for the three tables
table1_path = "concatenated_score_table.csv"
table2_path = "concatenated_score_table_papilio.csv"
table3_path = "concatenated_score_table_elegans.csv"

# Read in the data from each table as a pandas dataframe
df1 = pd.read_csv(table1_path)
df2 = pd.read_csv(table2_path)
df3 = pd.read_csv(table3_path)

# Extract the column of data from each dataframe
data1 = df1["SEQUENCE_SCORE"]
data2 = df2["SEQUENCE_SCORE"]
data3 = df3["SEQUENCE_SCORE"]

# Create a figure and axes object
fig, ax = plt.subplots()

# Set the title and y-axis label
ax.set_title("Guidance Scores for Arcellinida, H. papilio and H. elegans-specific genes")
ax.set_ylabel("Guidance Scores")

# Create the boxplots and add them to the axes
ax.boxplot([data1, data2, data3], labels=["Arcellinida", "H. papilio", "H. elegans"], showfliers=False)
ax.set_ylim([0, None])

# Show the plot
plt.show()

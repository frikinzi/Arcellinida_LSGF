library(sqldf)
library(dplyr)

AA_table = read.table(file = 'LSGF_Sequences/All_AA_Results.tsv', sep = '\t', header = TRUE)

NTD_table = read.table(file = 'LSGF_Sequences/All_NTD_Results.tsv', sep = '\t', header = TRUE)

res <- sqldf('SELECT query FROM NTD_table EXCEPT SELECT query FROM AA_table')

filtered_NTD <- dplyr::filter(NTD_table, query %in% res$query)

View(filtered_NTD)

write.table(filtered_NTD, "unique_BLAST_res.tsv", sep = '\t', row.names = FALSE)

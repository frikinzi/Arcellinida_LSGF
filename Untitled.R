library(tidyverse)

relativelength = read.table(file = 'hpap.csv', sep = ',', header = TRUE)
ggplot(relativelength, aes(color = OG, x=relative_length, y=as.factor(no_paralogs))) + geom_point() + theme(legend.position="none")

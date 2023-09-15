library(dplyr)
library(ggplot2)

setwd('/Users/brianamullins/Desktop')
#transcription sites
A <- read.table(file = 'tss_in_A_compartments.tsv', sep = '\t', header = TRUE)
B <- read.table(file = 'tss_in_B_compartments.tsv', sep = '\t', header = TRUE)

#ES is HK327ac peaks 
A_peaks <- read.table(file = 'peaks_ES_in_A_compartments.tsv', sep = '\t', header = TRUE)
B_peaks <- read.table(file = 'peaks_ES_in_B_compartments.tsv', sep = '\t', header = TRUE)

#HouseKeeping genes
A_House <- read.table(file = 'HK_genes_in_A_compartments.tsv', sep = '\t', header = TRUE)
B_House <- read.table(file = 'HK_genes_in_B_compartments.tsv', sep = '\t', header = TRUE)

sum(A$X107)
sum(B$X27)

sum(A_peaks$X76)
sum(B_peaks$X0)

sum(A_House$X9)
sum(B_House$X0)

CompartmentA <- c(24796, 25745, 3084)
CompartmentB <- c(11955, 3058, 252)
Genes <- c("Transcription_Start_Site", "Peaks", "House_Keeping_Genes")
genes_by_compartment <- data.frame(Genes, CompartmentA, CompartmentB)

p_A<-ggplot(data=genes_by_compartment, aes(x=Genes, y=CompartmentA)) +
  geom_bar(stat="identity", fill="steelblue") + ylim(0, 26000)
print(p_A + labs(title= "Density by Compartment A",
                 y="Density", x="Marker"))

p_B<-ggplot(data=genes_by_compartment, aes(x=Genes, y=CompartmentB)) +
  geom_bar(stat="identity", fill="steelblue") + ylim(0, 3000)
print(p_B + labs(title= "Denstiy by Compartment B",
                 y="Density", x="Marker"))

####### By Chromosome ########
pA<-ggplot(data=A, aes(x=X107, y=chr7)) +
  geom_bar(stat="identity", fill="steelblue") + xlim(0,2500) + coord_flip()
print(pA + labs(title= "Transcription Site Density by Chromosome (A)",
                 y="Chromosome", x="Density"))

pB<-ggplot(data=B, aes(x=X27, y=chr7)) +
  geom_bar(stat="identity", fill="steelblue") + xlim(0,2500) + coord_flip()
print(pB + labs(title= "Transcription Site Density by Chromsome (B)",
                y="Chromosome", x="Density"))

pA_peaks<-ggplot(data=A_peaks, aes(x=X76, y=chr7)) +
  geom_bar(stat="identity", fill="steelblue") + xlim(0,2200) + coord_flip()
print(pA_peaks + labs(title= "Histone Peak Density by Chromosome (A)",
                y="Chromosome", x="Density"))

pB_peaks<-ggplot(data=B_peaks, aes(x=X0, y=chr7)) +
  geom_bar(stat="identity", fill="steelblue") + xlim(0,2200) + coord_flip()
print(pB_peaks + labs(title= "Histone Peak Density by Chromosome (B)",
                      y="Chromosome", x="Density"))

pA_House<-ggplot(data=A_House, aes(x=X9, y=chr7)) +
  geom_bar(stat="identity", fill="steelblue") + xlim(0,300) + coord_flip()
print(pA_House + labs(title= "Housekeeping Gene Density by Chromosome (A)",
                      y="Chromosome", x="Density"))

pB_House<-ggplot(data=B_House, aes(x=X0, y=chr7)) +
  geom_bar(stat="identity", fill="steelblue") + xlim(0,300) + coord_flip()
print(pB_House + labs(title= "Housekeeping Gene Density by Chromosome (B)",
                      y="Chromosome", x="Density"))













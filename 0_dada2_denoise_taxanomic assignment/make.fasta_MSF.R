rm(list=ls())
library(dplyr)

#import taxa file
taxa <- as.data.frame(read.csv("C:/Users/jacqu/Desktop/BIOSSCOPE_SizeFrac_dada2_output/dada2 output/BIOSSCOPE_SizeFrac_taxa.csv"))

#subset SAR11, SAR202, and cyanobacteria ASVs
#Remove rows with chloroplast as the Order name
SAR11_taxa <- taxa %>% filter(Order =="SAR11 clade")
SAR202_taxa <- taxa %>% filter(Order =="SAR202 clade")
cyano_taxa <- taxa %>% filter(Phylum == "Cyanobacteria")

#start by giving our seq headers more manageable names (ASV_1, ASV_2...)
colnames <- c("ASV","Kingdom","Phylum","Class","Order","Family","Genus","Species")
colnames(sar11_taxa) <- colnames
colnames(sar202_taxa) <- colnames
colnames(cyano_taxa) <- colnames
sar11_seqs <- sar11_taxa$ASV
sar202_seqs <- sar202_taxa$ASV
pro_seqs <- cyano_taxa$ASV

sar11_headers <- vector(dim(sar11_taxa)[1], mode="character")
for (i in 1:dim(sar11_taxa)[1]) {
  sar11_headers[i] <- paste(">ASV", i, sep="_")
}

sar202_headers <- vector(dim(sar202_taxa)[1], mode="character")
for (i in 1:dim(sar202_taxa)[1]) {
  sar202_headers[i] <- paste(">ASV", i, sep="_")
}

pro_headers <- vector(dim(cyano_taxa)[1], mode="character")
for (i in 1:dim(cyano_taxa)[1]) {
  pro_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
sar11_fasta <- c(rbind(sar11_headers, sar11_seqs))
sar202_fasta <- c(rbind(sar202_headers, sar202_seqs))
pro_fasta <- c(rbind(pro_headers, pro_seqs))

write(sar11_fasta, "Baratheon_SAR11.fasta")
write(sar202_fasta, "Baratheon_SAR202.fasta")
write(pro_fasta, "Baratheon_Cyano.fasta")

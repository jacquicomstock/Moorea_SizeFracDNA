library(vegan)
library(tidyverse)

#remove eukaryotes, chloroplasts and mitochondria
setwd("C:/Users/jacqu/Desktop/Research/Projects/Moorea/MOOREA 2019/Baratheon/0_CLEANED.REDONE.analyses")
count_tab <- read.csv("MOOREA2019_Baratheon_seqtab-nochimtaxa.csv", header=TRUE,sep=",",row.names=1)
tax_tab <- read.csv("MOOREA2019_Baratheon_taxa.csv", header=TRUE,sep=",",row.names=1)
meta <- read.csv("MOOREA2019_Baratheon_metadata.csv", header=TRUE,sep=",",row.names=1)

#Remove rows with chloroplast as the Order name
ASV_nochl <- count_tab %>% filter(Order !="Chloroplast")
TAX_nochl <- tax_tab %>% filter(Order !="Chloroplast")
#Remove rows with mitochondria as the Family name
ASV_nochl.mit <- ASV_nochl %>% filter(Family !="Mitochondria")
TAX_nochl.mit <- TAX_nochl %>% filter(Family !="Mitochondria")
#Remove rows with Eukaryotes as the Kingdom name
ASV_nochl.mit.euk <- ASV_nochl.mit %>% filter(Kingdom !="Eukaryota")
TAX_nochl.mit.euk <- as.matrix(TAX_nochl.mit %>% filter(Kingdom !="Eukaryota"))

#write new files
write.csv(ASV_nochl.mit.euk, file="MOOREA2019_Baratheon_seqtab-nochimtaxa_NoChlMitEuk.csv", row.names=T)
write.csv(TAX_nochl.mit.euk, file="MOOREA2019_Baratheon_taxa_NoChlMitEuk.csv", row.names=T)

#generate rarefaction curve
ASV <- ASV_nochl.mit.euk[,1:(ncol(ASV_nochl.mit.euk)-7)]
tASV <- as.data.frame(t(ASV))
rarecurve(tASV, step = 1000, col = "blue", label=F, xlim=c(0,20000))

#sum reads in each sample to determine rarefaction depth
sum <- as.data.frame(colSums(ASV))

#rarefy samples to 19,000
OTU <- otu_table(ASV, taxa_are_rows = TRUE)
TAX <- tax_table(TAX_nochl.mit.euk)
phy <- phyloseq(OTU,TAX)

set.seed(8800)
rar <- rarefy_even_depth(phy, sample.size = 19000)
ASV_rar <- rar@otu_table
TAX_rar <- rar@tax_table

#calculate relative abundance
ASV_rar_prop <- apply(ASV_rar, 2, function(x) x/sum(x)*100)

#remove singletons (1/19000*100=0.0053%)
single <- rowSums(ASV_rar_prop[,1:ncol(ASV_rar_prop)]) > 0.0053
ASV_NoSingle <- as.data.frame(ASV_rar_prop[single,])

#Remove singletons from rarefied taxonomy file
ASVrow <- as.data.frame(rownames(ASV_NoSingle))
colnames(ASVrow) <- "ASV"
taxa_rar <- as.data.frame(TAX_rar)
taxa_rar$ASV <- rownames(taxa_rar)
taxa_NoSingle <- inner_join(taxa_rarA,ASVrow, by="ASV")

#create trimmed/rarefied metadata table
sample_rar <- as.data.frame(colnames(ASV_NoSingle))
colnames(sample_rar) <- "Sample"
sample_rar$Sample <- gsub('_','.', sample_rar$Sample)

meta$Sample <- rownames(meta)
meta$Sample <- gsub('_','.', meta$Sample)
meta_all <- full_join(sample_rar,meta, by= "Sample")
meta_rar <- inner_join(sample_rar, meta, by="Sample")

#write new ASV and taxa file
write.csv(taxa_NoSingle, file="Baratheon_taxa_rar19000_NoSingle.csv")
write.csv(ASV_NoSingle, file="Baratheon_seqtab-nochimtaxa_rar19000_NoSingle.csv")
write.csv(meta_rar,file="Baratheon_meta_rar19000.csv")

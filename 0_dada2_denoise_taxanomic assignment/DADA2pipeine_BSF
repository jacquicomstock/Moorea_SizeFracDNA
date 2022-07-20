library("dada2")
path <- "/home/jcomstock/BIOS_Frac"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    truncLen=c(230,160), 
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE, matchIDs=TRUE)
}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out
saveRDS(errF, "errF.rds")
saveRDS(errR, "errR.rds")
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaFs, "dadaFs_N.rds")
saveRDS(dadaRs, "dadaRs_N.rds")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
saveRDS(mergers, "dadaFs_N.rds")
seqtab <- makeSequenceTable(mergers[names(mergers) != "Mock"])
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
saveRDS(seqtab.nochim, "dadaFs_Nseqtab.nochim.rds")
t.seqtab.nochim <- t(seqtab.nochim)
asv_seqs <- rownames(t.seqtab.nochim)
asv_headers <- vector(dim(t.seqtab.nochim)[1], mode="character")
for (i in 1:dim(t.seqtab.nochim)[1]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "BIOSSCOPE_SizeFrac_ASVs.fasta")
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_wSpecies_train_set.fa.gz")
unname(head(taxa))
saveRDS(taxa, "taxa.rds")
write.table(cbind(t(seqtab.nochim) , taxa), "BIOSSCOPE_SizeFrac_seqtab-nochimtaxa.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(taxa,"BIOSSCOPE_SizeFrac_taxa.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
q()

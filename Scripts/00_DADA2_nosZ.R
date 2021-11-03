
library(pacman)
p_load(dada2, tidyverse,phyloseq, Biostrings, here, DECIPHER, phangorn, ape, phylogram)
path <- "~/Documents/SaltikovLab/Sequencing/MAR_All_Sites/nosZ/"
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2]) 
plotQualityProfile(fnRs[1:2])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
  
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,200),
                      maxN=0, maxEE=c(2,2), trimLeft = c(23,21), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "Objects/nosZ_seqtabnochim.rds")




#make taxonomy datbase:
#run  accession_convert.sh in bash to get nosZ_db.fa which was made by pulling seed sequences from fungene
#zsh ~/Documents/SaltikovLab/Sequencing/MAR_All_Sites/nosZ/accession_convert.sh

taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/SaltikovLab/Sequencing/MAR_All_Sites/nosZ/nosZ_db.fa", multithread=TRUE)

saveRDS(taxa, here("Objects/nosZ_taxa.rds"))

seqtab.nochim<-readRDS(here("Objects/nosZ_seqtabnochim.rds"))
taxa<- readRDS(here("Objects/nosZ_taxa.rds"))

#make metadata file
Location = c("KTR","KTR","KTR","KTR","KTYA","KTYA","KTYA","KTYA")
Timing = c("Before","After","Before","After","Before","After","Before","After")
Treatment= c("NS","NS","WC","WC","NS","NS","WC","WC")

metadata<- tibble(Sample = rownames(seqtab.nochim),
                  Location = factor(Location),
                  Treatment= factor(Treatment),
                  Timing = factor(Timing, levels = c("Before", "After"), ordered = T ))%>%
  column_to_rownames(.,var="Sample")%>%
  mutate(., Type = paste0(Treatment, ": ", Timing))


#make phyloseq object
ps<- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
              sample_data(metadata),
              tax_table(taxa))

genus_ps<- tax_glom(ps, taxrank = "Genus", NArm = F)


#alignment then phylogentic tree
#get ref seq

referenceseq<-"~/Documents/SaltikovLab/Sequencing/MAR_All_Sites/nosZ/refsequences.fasta"
rseq<-readDNAStringSet(referenceseq)

#get sequences from glommed phyloseq
dna <- Biostrings::DNAStringSet(taxa_names(genus_ps))

names(dna) <- taxa_names(genus_ps)
#merge to get refseq populated
genus_ps <- merge_phyloseq(genus_ps, dna)
#change ASV to numbers
taxa_names(genus_ps) <- paste0("ASV", seq(ntaxa(genus_ps)))
#pull out the seq tab and taxa tab for making a new phyloseq
seqtab.nochim_genus<- otu_table(genus_ps)
taxa_genus<-tax_table(genus_ps)
#reset the names to ASV#
names(dna)<-taxa_names(genus_ps)

#make alignment
alignment <- AlignSeqs(c(dna, rseq))
  #saveRDS(alignment, here("Objects/nosZ_alignment.rds"))

phang.align <- phyDat(as(alignment, "matrix"), type="DNA") 
#saveRDS(phang.align,here("Objects/nosZ_phang_align.rds"))
dm <- dist.ml(phang.align)
#aveRDS(dm, here("Objects/dm_nosZ.rds"))
treeNJ <- NJ(dm) # Note, tip order != sequence order
#saveRDS(treeNJ, here("Objects/treeNJ_nosZ.rds"))
fit = pml(treeNJ, data=phang.align)
#saveRDS(fit, here("Objects/fit_nosZ.rds"))
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

saveRDS(fitGTR, "Objects/nosz_phytree.RDS")

genus_ps<-phyloseq(otu_table(seqtab.nochim_genus),
                   tax_table(taxa_genus),
                   sample_data(metadata),
                   phy_tree(fitGTR$tree))

saveRDS(genus_ps, here("Objects/nosZ_phylo.rds"))



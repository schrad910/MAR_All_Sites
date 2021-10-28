library(here)
library(tidyverse)
library(dada2)
#loading the files from DADA2 pipelines in the Run1(2,3) folders
table1<-readRDS(here("Objects/DADA2/seqtabnochim_run1.rds"))
table2<-readRDS(here("Objects/DADA2/seqtabnochim_run2.rds"))
table3<-readRDS(here("Objects/DADA2/seqtabnochim_run3.rds"))

#RE-MAKING TAXA FROM THE MERGED TABLES
#taxa1<-readRDS(here("Objects/DADA2/taxa_run1.rds"))
#taxa2<-readRDS(here("Objects/DADA2/taxa_run2.rds"))
#taxa3<-readRDS(here("Objects/DADA2/taxa_run3.rds"))

#combining taxa table and removing duplicates
#taxa<- rbind(taxa1,taxa2,taxa3)

#saveRDS(taxa, here("Objects/taxa.rds"))

#combining sequence table and saving
#table<-mergeSequenceTables(table1,table2,table3)

#saveRDS(table, here("Objects/table.rds"))

table<- readRDS(here("Objects/table.rds"))
## try 1: Using the big data tutorial and using the combined table  make a taxa table
#update:try 1 =sucessful
#taxa <- assignTaxonomy(table, "~/Documents/SaltikovLab/Sequencing/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)

#saveRDS(taxa, here("Objects/taxa.rds"))
taxa<- readRDS(here("Objects/taxa.rds"))

library(DECIPHER)
library(phangorn)
seqs <- getSequences(table)
names(seqs) <- seqs 

alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
saveRDS(alignment, here("Objects/alignment.rds"))
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
saveRDS(phangAlign,here("Objects/phangAlign.rds"))
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
saveRDS(fitGTR, here("Objects/fitGTR.rds"))
#bring in metadata file & set ojects as factors
library(phyloseq)
library(Biostrings)
metadata <-read_tsv(here("Environment/metadata.txt"), 
              col_types = cols(Timing=col_factor(levels = c("Before", "After")),
                      Treatment=col_factor(levels = c("NS","MX","WC","BC")),
                      Location=col_factor(levels=c("KTR","KTYA","HSP"))))%>%
  column_to_rownames(var="SampleID")%>%
  sample_data()
  #make phyloseq object

ps<- phyloseq(otu_table(table, taxa_are_rows = FALSE),
              sample_data(metadata),
              tax_table(taxa))
#renamimg to ASV1, ASV2, etc. 
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

aligngenus<- AlignSeqs(refseq(genus), anchor=NA)
gphangAlign <- phyDat(as(aligngenus, "matrix"), type="DNA")
dm <- dist.ml(gphangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=gphangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

saveRDS(ps, file = here("Objects/phlyoseq.rds"))


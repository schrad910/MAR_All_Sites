library("phyloseq")
library("DESeq2")

#combining lines with same genus
phylo<-tax_glom(phylo,taxrank = "Genus" )
#removing samples that have less than 3900 reads *one replicate has this
#phylo<-prune_samples(sample_sums(phylo)>3900,phylo)

trmt_deseq= phyloseq_to_deseq2(phylo, ~Timing+Treated)
#calculate geometric means before due to the large amounts of zeros in data
gm_mean= function(x, na.rm=TRUE) {
  exp(sum(log(x[x>0]), na.rm = na.rm) / length(x))
}
geomeans= apply(counts(trmt_deseq), 1, gm_mean)
trmt_deseq= estimateSizeFactors(trmt_deseq, geoMeans= geomeans)
#using Wald test because its commonly used
trmt_deseq = DESeq(trmt_deseq,test = "Wald", fitType = "parametric")
res= results(trmt_deseq, cooksCutoff = FALSE)
alpha=0.01
sigtab = res[which(res$padj< alpha),]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phylo)[rownames(sigtab),],"matrix"))
sigtab<- sigtab[order(sigtab$log2FoldChange),]


#joinablesigtab <- rownames_to_column(sigtab, var="featureid")

theme_set(theme_bw())
#filter out NA genuses
sigtabgen = subset(sigtab, !is.na(Genus) & Genus != 'uncultured')


normalized_relabun<-counts(trmt_deseq, normalized=TRUE)%>%
  as_tibble(rownames="featureID")%>%
  mutate_at(.,vars(-featureID), .funs = ~./sum(.)) 

#join to get taxonomy

norm_tax_rel_abun<-inner_join(taxonomy, normalized_relabun, by= c("Feature.ID"="featureID"))%>%
  select(-Feature.ID)

  
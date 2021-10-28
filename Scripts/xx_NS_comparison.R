library(phyloseq)
library(tidyverse)
library(ggplot2)
library(DESeq2)
#pulling in the phyloseq object that was glomed for genus because most species are "uncultured"
phylo <- readRDS("~/Documents/SaltikovLab/Sequencing/MAR_All_Sites/Objects/genera.rds")
sample_data(phylo)$Location <-factor(get_variable(phylo, "Location"), c("HSP", "KTR", "KTYA"), ordered = FALSE )

sample_data(phylo)$Timing <-factor(get_variable(phylo, "Timing"), c("Before", "After"), levels = c("Before", "After"), ordered = FALSE)
#subsetting phyloseq object to only be Native soil
NS_phylo <- subset_samples(phylo, Treatment=="NS")

KTYA_NS <- subset_samples(NS_phylo, Location == "KTYA" )
KTR_NS <- subset_samples(NS_phylo, Location == "KTR" )
HSP_NS <- subset_samples(NS_phylo, Location == "HSP" )


ds<- phyloseq_to_deseq2(NS_phylo, ~Location +  Timing)

#collapse by replicates which are mainly close to each other in NDMS 



#remove genes that have zero counts for all samples
ds <- ds[ rowSums( counts(ds) ) > 0 , ]
collpse_dds<-DESeq(ds)
col_res <- results(collpse_dds)
col_res <- col_res[order(col_res$padj, na.last=NA), ]
alpha = 0.01
col_sigtab = col_res[(col_res$padj < alpha), ]
col_sigtab = cbind(as(col_sigtab, "data.frame"), as(tax_table(NS_phylo)[rownames(col_sigtab), ], "matrix"))

col_posigtab = col_sigtab[col_sigtab[, "log2FoldChange"] > 0, ]

#Compare when not collapsed
ds<- ds[rowSums(counts(ds))>0,]
# every gene has a zero when not collapsed so need to replace 0 with 1
gm_mean<-function(x, na.rm=TRUE) {
  exp(sum(log(x[x>0]), na.rm = na.rm) / length(x))
}
geomeans<-apply(counts(ds), 1, gm_mean)
ds<-estimateSizeFactors(ds, geoMeans=geomeans)
ds<- DESeq(ds)
res<- results(ds)%>%
  .[order(res$padj, na.last = NA),]
sigtab<- res[(res$padj <alpha),]
sigtab <-cbind(as(sigtab,"data.frame"), as(tax_table(NS_phylo)[rownames(sigtab),],"matrix"))

collapse<- col_sigtab%>%
  select(Genus, log2FoldChange)
non_collapse <- sigtab %>%
  select(Genus, log2FoldChange)
comparison<- full_join(non_collapse, collapse, by= "Genus") %>%
  subset(!is.na(Genus) | Genus != "uncultured")

## looking at classes

#class<-readRDS("~/Documents/SaltikovLab/Sequencing/MAR_All_Sites/Objects/classes.rds")
NS_class<-  subset_samples(class, Treatment=="NS")

class_ds<- phyloseq_to_deseq2(NS_class, ~Location + Timing)



#remove genes that have zero counts for all samples
class_ds <- class_ds[ rowSums( counts(class_ds) ) > 0 , ]
class_ds<-DESeq(class_ds)
class_res <- results(class_ds)
class_res <- class_res[order(class_res$padj, na.last=NA), ]
class_sigtab = class_res[(class_res$padj < alpha), ]
class_sigtab = cbind(as(class_sigtab, "data.frame"), as(tax_table(NS_class)[rownames(class_sigtab), ], "matrix"))
class_sigtab
class_posigtab = class_sigtab[class_sigtab[, "log2FoldChange"] > 0, ]


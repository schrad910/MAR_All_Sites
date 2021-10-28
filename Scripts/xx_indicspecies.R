library(pacman)
set.seed(3618)
p_load(tidyverse, ggpubr, indicspecies, here)

#this code is looking at indicator species of before and after infiltratio
phylo<-readRDS(here("Objects/genera.rds"))
native_soil<-subset_samples(phylo, Treatment=="NS")
table<-as(otu_table(native_soil), "matrix")
taxa<-as(tax_table(native_soil)), "matrix")%>%
  as_tibble(rownames = "Label")
metadata<-as(sample_data(native_soil),"matrix")%>%
  as.tibble(rownames=NA)


timing<-metadata$Timing
location<-metadata$Location

#run indicspecies  999 permutations, func= point biserial correlation coefficient to look at timing
inv<-multipatt(table, timing, func = "r.g", restcomb=location, control = how(nperm=999))
summary(inv)

#run indicspecies to look at location

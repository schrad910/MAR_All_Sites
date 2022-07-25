library(pacman)
set.seed(3618)
p_load(tidyverse,phyloseq, ANCOMBC, here, knitr,ggpubr, vegan, ggvegan)
phylo<- readRDS(here("Objects/genera.rds"))
phylo<-subset_samples(phylo,Treatment %in% c("NS", "WC"))
table<-as(otu_table(phylo), "matrix")
ndms<-metaMDS(table, distance = "bray", try=40, autotransform = F)
metadata <- readRDS(here("Objects/metadata.rds"))%>%
  subset(., Treatment %in% c("NS", "WC"))
metadata$Location<- factor(metadata$Location, levels = c("HSP", "KTR", "KTYA"), ordered = T)
metadata<-mutate(metadata, replicate_group= paste0(Location," ",Treatment," ",Timing))
f_ndms<- fortify(ndms)%>%
  subset(., Score=="sites")%>%
  left_join(., metadata, by=c("Label"="SampleID"))
anova<- adonis2(table~replicate_group, data = metadata, 
                permutations = 10000 )
f_ndms$Location<-factor(f_ndms$Location, levels = c("HSP", "KTR", "KTYA"), ordered = T)
ggplot(subset(f_ndms))+
  geom_point(aes(x=NMDS1,y=NMDS2, col=replicate_group, shape=Location))+
  theme_pubr(legend = "right")+
  #scale_color_manual(values=palette)+
  #scale_fill_manual(values=palette)+
  scale_shape_manual(values = c(20,4, 25))+ 
  geom_abline(intercept = 0,slope = 0,linetype="dashed", size=0.8)+ 
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.8)

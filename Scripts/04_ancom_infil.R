#this code does ancom for TIMNG 
library(pacman)
set.seed(3618)
p_load(tidyverse,phyloseq, ANCOMBC, here, knitr,ggpubr, vegan)
phylo<- readRDS(here("Objects/genera.rds"))
native_soil<-subset_samples(phylo, Treatment=="NS")
#filtering out taxa not seen more than 10 times in 25% of the samples to reduce low abundance 
native_soil<-filter_taxa(native_soil, function(x) sum(x >10) > (0.25*length(x)), TRUE)
taxa<-as(tax_table(native_soil), "matrix")%>%
  as_tibble(rownames = "Label")
table<-as(otu_table(native_soil), "matrix")
#doing ancom-bias correction for infiltration, form= how absoulte values depend on loaction ad timing
#taxa w/ .90% zeros across sampels are cut, no rarefication, group is timing, 
#
out = ancombc(phyloseq= native_soil, formula = "Timing"  , p_adj_method = "holm", zero_cut=.9,
              lib_cut = 0, struc_zero=FALSE, neg_lb = FALSE,  global=FALSE, alpha=0.05) 

#make a data frame of results
out_df <- data.frame(
  Species = row.names(out$res$beta),
  coef = unlist(out$res$beta),
  se = unlist(out$res$se),
  W = unlist(out$res$W),
  p_val = unlist(out$res$p_val),
  q_val = unlist(out$res$q_val),
  diff_abn = unlist(out$res$diff_abn))
#filter only  significant genera and change from natural log to 2-log
out_df<- subset(out_df, diff_abn==TRUE)%>%
  mutate(., log2= coef * log2(exp(1)))%>%
  mutate(., selog2= se *log2(exp(1)))%>%
  left_join(., taxa, by= c("Species" = "Label"))
#change "weird" uncultured genus names to family name
out_df<- mutate(out_df, Genus= if_else(grepl('[0-9]', Genus),paste0(Family," ", Genus), Genus))%>%
  mutate(., Genus= if_else(Genus=="hgcI clade", paste0(Family," ", Genus) , Genus))
#saveRDS(out_df, here("Objects/out_df.rds"))
palette<-colorRampPalette(readRDS(here("Objects/palette.RDS")))(9)

ggplot(data=out_df)+
  geom_col(aes(x=log2, y=reorder(Genus,-log2), fill=Phylum))+
  geom_errorbarh(aes(y=Genus, xmin= log2-selog2,xmax=log2+selog2), height=.2)+
  scale_fill_manual(values = palette)+
  theme_pubr(legend="right")+
  labs(x= "Log 2-Fold Abundance Difference", y= "Genus")+
  theme(legend.text = element_text(face = "italic"), axis.text.y  = element_text(face = "italic"))
ggsave(here("Figures/Figures/Figure04.tiff"))  
#checking with glommed phylums

#class_ns<- tax_glom(native_soil, taxrank = "Phylum", NArm = F)

#class_out<-ancombc(phyloseq= class_ns, formula = "Timing", p_adj_method = "holm", zero_cut=.9,
 #                   lib_cut = 0, struc_zero=FALSE, neg_lb = FALSE,  global=FALSE, alpha=0.05) 
#class_df <- data.frame(
  #Species = row.names(class_out$res$beta),
#  coef = unlist(class_out$res$beta),
#  se = unlist(class_out$res$se),
#  W = unlist(class_out$res$W),
#  p_val = unlist(class_out$res$p_val),
#  q_val = unlist(class_out$res$q_val),
 # diff_abn = unlist(class_out$res$diff_abn))
#class_df<-subset(class_df, diff_abn==TRUE)%>%
 # left_join(., taxa, by= c("Species" = "Label"))%>%
#  select(1:12)
#class_df<-mutate(class_df, log2= coef* log2(exp(1)))%>%
 # mutate(., selog2= se * log2(exp(1)))

#ggplot(data=class_df)+
 # geom_col(aes(x=log2, y=reorder(Phylum,-log2), fill=Phylum))+
  #geom_errorbarh(aes(y=Phylum, xmin= log2-selog2,xmax= log2+selog2), height=.2)+
  #scale_fill_manual(values = palette)+
  #theme_pubr(legend="right")+
  #coord_fixed()+
  #labs(y="Class", x= "Log 2-fold Change")
#ggsave(here("Figures/ancom_timing_class.png"))         



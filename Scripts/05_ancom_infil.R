#this code does ancom for TIMNG 
library(pacman)
set.seed(3618)
p_load(tidyverse,phyloseq, ANCOMBC, here, knitr,ggpubr, vegan)
phylo<- readRDS(here("Objects/genera.rds"))
native_soil<-subset_samples(phylo, Treatment=="NS")
#native_soil<-tax_glom(native_soil, taxrank = "Order")
#filtering out taxa not seen more than 10 times in 25% of the samples to reduce low abundance 
native_soil<-filter_taxa(native_soil, function(x) sum(x >100) > (0.25*length(x)), TRUE)

taxa<-as(tax_table(native_soil), "matrix")%>%
  as_tibble(rownames = "Label")
table<-as(otu_table(native_soil), "matrix")
meta<-as(sample_data(native_soil),"matrix")

ancom_infil<- function(phyloseq) {

#doing ancom-bias correction for infiltration, form= how absoulte values depend on loaction ad timing
#taxa w/ .90% zeros across sampels are cut, no rarefication, group is timing, 
#
  out = ancombc(phyloseq= phyloseq, formula = "Timing"  , p_adj_method = "holm", zero_cut=.9,
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
}

palette<-colorRampPalette(p)(13)
p<-c("#89C0DB", "#593F28" ,"#739D6E" ,"#D0C1B0" ,"#80654C","#70a08c")
  # ggplot(data=out_df)+
  #   geom_col(aes(x=log2, y=reorder(Genus,-log2), fill=Phylum))+
  #   geom_errorbarh(aes(y=Genus, xmin= log2-selog2,xmax=log2+selog2), height=.2)+
  #   scale_fill_manual(values = palette)+
  #   theme_pubr(legend="right")+
  #   labs(x= "Log 2-Fold Abundance Difference", y= "Genus")+
  #   theme(legend.text = element_text(face = "italic"), axis.text.y  = element_text(face = "italic"))
  # 


native_soil_hsp<-subset_samples(native_soil, Location=="HSP")
native_soil_ktr<-subset_samples(native_soil, Location=="KTR")
native_soil_ktya<-subset_samples(native_soil, Location=="KTYA")
hsp<-ancom_infil(native_soil_hsp)%>%
  mutate(.,location="HSP")
ktr<-ancom_infil(native_soil_ktr)%>%
  mutate(.,location="KTR")
ktya<-ancom_infil(native_soil_ktya)%>%
  mutate(.,location="KTYA")
#checking to make sure results don't change with only one replicate being used
# native_soil_hsp1<-subset_samples(native_soil, Location=="HSP" & Replicate==1)
# native_soil_ktr1<-subset_samples(native_soil, Location=="KTR"& Replicate==1)
# native_soil_ktya1<-subset_samples(native_soil, Location=="KTYA"& Replicate==1)
# hsp1<-ancom_infil(native_soil_hsp1)%>%
#   mutate(.,location="HSP")
# ktr1<-ancom_infil(native_soil_ktr1)%>%
#   mutate(.,location="KTR")
# ktya1<-ancom_infil(native_soil_ktya1)%>%
#   mutate(.,location="KTYA")
# identical(hsp,hsp1)
# identical(ktr, ktr1)
# identical(ktya,ktya1)
combined<-rbind(hsp,ktr,ktya)


palette<-c("#89c0db","#2c7dbc","#003790","#cbe2ff","#8092bb",
           "#615e90","#e2d3c1",
           "#b29f8a","#856d56","#593f28","#97d0a4","#447958","#004429")
#saveRDS(palette,file=here("Objects/palette.RDS"))
ggplot(data=combined)+
  geom_col(aes(x=log2, y=reorder(Genus,-log2), fill=Phylum))+
  geom_errorbarh(aes(y=Genus, xmin= log2-selog2,xmax=log2+selog2), height=.2)+
  geom_vline(aes(xintercept=0),linetype=2, color="grey")+
  facet_grid(cols = vars(location))+
  scale_fill_manual(values = palette)+
  theme_pubr(legend="right")+
  labs(x= "Log 2-Fold Abundance Difference", y= "Genus", caption= "Native soil (no PRB)")+
  theme(legend.text = element_text(face = "italic"), 
        axis.text.y  = element_text(face = "italic"),
        axis.title.y = element_text(vjust=1.02,hjust=1, face="bold", angle=0))


ggsave(here("round3/Figure05.tiff"),height = 9, unit="in")  

inner_join(hsp, ktya, by="Species")%>%
  select(Genus.x,log2.x,log2.y)
#8 commonalities between HSP and KTYA!  17 showed up with both but 9 were opposites

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



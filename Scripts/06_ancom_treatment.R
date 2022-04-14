library(pacman)
set.seed(3618)
p_load(tidyverse,phyloseq, ANCOMBC, here, knitr,ggpubr, vegan, ggvegan)
phylo<- readRDS(here("Objects/genera.rds"))
#looking at just WC vs NS changes
metadata <- readRDS(here("Objects/metadata.rds"))%>%
  subset(., Timing == "After" & Treatment%in%  c("WC","NS"))
after<- subset_samples(phylo, Timing == "After" & Treatment %in%  c("WC","NS"))
#filter outtaxa not seen more than 10 times in 25% of the samples to reduce low abundance 
after<- filter_taxa(after, function(x) sum(x > 52) > (0.25*length(x)), TRUE)
taxa<-as(tax_table(after), "matrix")%>%
  as_tibble(rownames = "Label")
table<-as(otu_table(after), "matrix")



out = ancombc(phyloseq= after, formula = "Treatment", p_adj_method = "holm", zero_cut=.9,
               struc_zero=FALSE,   alpha=0.05) 
out_df <- data.frame(
  Species = row.names(out$res$beta),
  coef = unlist(out$res$beta),
  se = unlist(out$res$se),
  W = unlist(out$res$W),
  p_val = unlist(out$res$p_val),
  q_val = unlist(out$res$q_val),
  diff_abn = unlist(out$res$diff_abn))

out_df<- subset(out_df, diff_abn==TRUE)%>%
  mutate(., log2= coef * log2(exp(1)))%>%
  mutate(., selog2= se *log2(exp(1)))%>%
  left_join(., taxa, by= c("Species" = "Label"))

out_df<- mutate(out_df, Genus= if_else(grepl('[0-9]', Genus),paste0(Family," ", Genus), Genus))%>%
  mutate(., Genus= if_else(Genus=="hgcI clade", paste0(Family," ", Genus), Genus))



palette<-colorRampPalette(readRDS(here("Objects/palette.RDS")))(8)
palette1<- c("#7D9FAE" ,"#593F28" ,  "#89C0DB" ,"#80654C",  "#476129" ,  "#655F54" ,
              "#A8937E" , "#8AA67E" ,"#B8B89F" ,"#6C855C"  )
ggplot(data=out_df, aes(x=log2, y=reorder(Genus, -log2), fill=Phylum))+
  geom_col( position = "dodge")+
  geom_errorbar(aes( xmin= log2-selog2,xmax=log2+selog2), width=.2, position = position_dodge(1))+
  scale_fill_manual(values = palette1)+
  theme_pubr(legend="right")+
  labs(x= "Log 2- Fold Abundance Difference", y= "Genus")+
  coord_fixed()+  
  theme(legend.text = element_text(face = "italic"), axis.text.y  = element_text(face = "italic"))

  ggsave(here("Figures/Round2/Figure06.tiff"), width = 7, height = 7, units = "in")



#looking at WC after compared to before
#WC_metadata <- readRDS("~/Documents/SaltikovLab/Sequencing/MAR_All_Sites/Objects/metadata.rds")%>%
#  subset(.,Treatment=="WC")
#WC_after<- subset_samples(phylo, Treatment =="WC")
#filter outtaxa not seen more than 10 times in 25% of the samples to reduce low abundance 
#WC_after<- filter_taxa(WC_after, function(x) sum(x > 10) > (0.25*length(x)), TRUE)
#WC_taxa<-as(tax_table(WC_after), "matrix")%>%
 # as_tibble(rownames = "Label")
#WC_table<-as(otu_table(WC_after), "matrix")



#WC_out = ancombc(phyloseq= WC_after, formula = "Timing", p_adj_method = "holm", zero_cut=.9,
 #             struc_zero=FALSE,   alpha=0.05) 
#WC_out_df <- data.frame(
#  Species = row.names(WC_out$res$beta),
# coef = unlist(WC_out$res$beta),
 # se = unlist(WC_out$res$se),
#  W = unlist(WC_out$res$W),
 # p_val = unlist(WC_out$res$p_val),
 # q_val = unlist(WC_out$res$q_val),
# diff_abn = unlist(WC_out$res$diff_abn))

#WC_out_df<- subset(WC_out_df, diff_abn==TRUE)%>%
 # mutate(., log2= coef * log2(exp(1)))%>%
  #mutate(., selog2= se *log2(exp(1)))%>%
  #left_join(., WC_taxa, by= c("Species" = "Label"))

#WC_out_df<- mutate(WC_out_df, Genus= if_else(grepl('[0-9]', Genus),Family, Genus))
#WC_out_df<-mutate(WC_out_df, Genus= if_else(Genus=="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobium", Genus))


#b<-ggplot(data=WC_out_df, aes(x=log2, y=reorder(Genus, -log2), fill=Phylum))+
 # geom_col( position = "dodge")+
  #geom_errorbar(aes( xmin= log2-selog2,xmax=log2+selog2), width=.2, position = position_dodge(1))+
  #scale_fill_manual(values = c(palette))+
  #theme_pubr(legend="right")+
  #labs(x= "Log 2- Fold Abundance Difference", y= "Genus")+
  #coord_fixed()

ns_out_df<-readRDS(here("Objects/out_df.rds"))
join<-inner_join(ns_out_df, out_df, by="Species")%>%
  select(2,15,16)

join<- rename(join, "Genus"="Genus.x", "Infiltration Changes"= "coef.x",  "Treatment Changes"="coef.y")%>%
  select(Genus, 'Infiltration Changes', 'Treatment Changes')



#saving code for when there is 2+ factros


#c_out_df <- data.frame(res$beta * res$diff_abn, check.names = F)%>%
  # mutate(., total= TreatmentMX+ TreatmentWC + TreatmentBC)%>%
 # subset(., TreatmentMX!= 0 & TreatmentWC!= 0 ) %>%
  #select(-total)%>%
  #rename(., "WC"= "TreatmentWC", "MX"="TreatmentMX")%>%
  #rownames_to_column(var="ASV")

#c_out_df<- pivot_longer(c_out_df, cols = c("WC","MX"), names_to = "Treatment", values_to = "coef"  )

#c_SE<- data.frame(res$se * res$diff_abn, check.names = F)%>%
 # mutate(., total= TreatmentMX+ TreatmentWC)%>%
  #subset(., total!= 0) %>%
 # select(-total)%>%
 # rename(., "WC"= "TreatmentWC", "MX"="TreatmentMX")%>%
  #rownames_to_column(var="ASV")


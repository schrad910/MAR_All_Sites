library("phyloseq")
library("DESeq2")
library("tidyverse")
library("here")

#glommed to genus and saved.... may have to go up in ranks
genus<- readRDS(here("Objects/genera.rds"))
#saveRDS(genus, "~/Documents/SaltikovLab/Sequencing/MAR_All_Sites/Objects/genera.rds")

#first looking at NS samples to see the infiltration effect

infiltration_ps<-subset_samples(genus, Treatment == "NS")
infil_ds<-phyloseq_to_deseq2(infiltration_ps, ~ Timing )

#making DESEQ object
infil_ds<-DESeq(infil_ds)

#pulling table of normalized counts
normalized_infil<-counts(infil_ds,normalized=TRUE)%>%
  as_tibble(rownames = "ASV")%>%
 mutate_at(.,vars(-ASV),.funs = ~./sum(.))
alpha=.01
#trying agin w/ relevl
infil_ds$Location<- relevel(infil_ds$Location, "KTYA")
infil_ds<-DESeq(infil_ds)
res1<- results(infil_ds)
res1<-results(infil_ds, contrast = c("Timing", "After","Before"))
sigres1<-res1[which(res1$padj<alpha),]
sigres1<- cbind(as(sigres1,"data.frame"),as(tax_table(genus)[rownames(sigres1),],"matrix"))
sigres1<- sigres1[order(sigres1$log2FoldChange),]

ggplot(sigres1, aes(x=log2FoldChange, y=reorder(Genus, -log2FoldChange), fill=Phylum))+
  geom_col()+
  geom_vline(xintercept = 0,aes(size=2))

sigres2<-res2[which(res2$padj<alpha),]
sigres2<- cbind(as(sigres2,"data.frame"),as(tax_table(genus)[rownames(sigres2),],"matrix"))
sigres2<- sigres2[order(sigres2$log2FoldChange),]
res3<-results(infil_ds, contrast =list(c("Timing_After_vs_Before","LocationHSP.TimingAfter")))
sigres3<-res3[which(res3$padj<alpha),]
sigres3<- cbind(as(sigres3,"data.frame"),as(tax_table(genus)[rownames(sigres3),],"matrix"))
sigres3<- sigres3[order(sigres3$log2FoldChange),]
palette<- readRDS(here("Objects/palette.RDS"))
#effect of infiltration for KTR (main effect)
res_infil_ktr<-results(infil_ds, contrast = c("Timing", "After","Before"))
sig_res_infil_ktr<- res_infil_ktr[which(res_infil_ktr$padj<alpha),]
sig_res_infil_ktr<- cbind(as(sig_res_infil_ktr,"data.frame"),as(tax_table(genus)[rownames(sig_res_infil_ktr),],"matrix"))
sig_res_infil_ktr<- sig_res_infil_ktr[order(sig_res_infil_ktr$log2FoldChange),]

#effect of infiltration for KTYA
#res_infil_ktya<- results(infil_ds, contrast = list(c("Timing_After_vs_Before","LocationKTYA.TimingAfter")))
#sig_res_infil_ktya<- res_infil_ktya[which(res_infil_ktya$padj<alpha),]
#sig_res_infil_ktya<- cbind(as(sig_res_infil_ktya,"data.frame"),as(tax_table(genus)[rownames(sig_res_infil_ktya),],"matrix"))
#sig_res_infil_ktya<- sig_res_infil_ktya[order(sig_res_infil_ktya$log2FoldChange),]

#effect of infiltration for HSP
#res_infil_hsp<- results(infil_ds, contrast = list(c("Timing_After_vs_Before","LocationHSP.TimingAfter")))
#sig_res_infil_hsp<- cbind(as(sig_res_infil_hsp,"data.frame"),as(tax_table(genus)[rownames(sig_res_infil_hsp),],"matrix"))
#sig_res_infil_hsp<- sig_res_infil_hsp[order(sig_res_infil_hsp$log2FoldChange),]


#comparisons between results:
#comparison from affect from 3 locations
loc_comp <- tibble(ASV=row.names(sigres1),LFC= sigres1$log2FoldChange,Phylum=sigres1$Phylum, Genus=sigres1$Genus, loc="KTYA")%>%
  rbind(tibble(ASV=row.names(sigres2),LFC= sigres2$log2FoldChange,Phylum=sigres2$Phylum, Genus=sigres2$Genus, loc="KTR"))%>%
  rbind(tibble(ASV=row.names(sigres3),LFC= sigres3$log2FoldChange,Phylum=sigres3$Phylum, Genus=sigres3$Genus, loc="HSP"))
  
#comparison between (KTYA+HSP nd just KTYA:b4 and after)
#form_comp<-tibble(ASV=row.names(sig_res_infil_ktya),LFC= sig_res_infil_ktya$log2FoldChange,Genus=sig_res_infil_ktya$Genus, base="KTR")%>%
 # rbind(tibble(ASV=row.names(sigres1),LFC= sigres1$log2FoldChange,Genus=sigres1$Genus, base="KYTA"))
#form_dif_comp<-c("ASV1011", "ASV16707","ASV1880", "ASV3000", "ASV4095","ASV4448", "ASV4458", "ASV7084","ASV824","ASV9287")
#form_unq_comp<-c("ASV128","ASV4656","ASV9541")

  
#making change fold difference chart split out by location
ggplot(loc_comp, aes(x=LFC, y=Genus, fill=Phylum))+
  geom_col()+
  facet_grid(cols= vars(loc))+
  geom_vline(xintercept = 0,aes(size=2))

#filtering to only see genus affected at all locations

f_loc_com<- pivot_wider(data = loc_comp,id_cols = c("ASV", "Phylum","Genus"), names_from = loc, values_from = LFC)%>%
  subset(!is.na(KTR) & !is.na(HSP))%>%
  pivot_longer(cols = c(KTYA,KTR,HSP), names_to= "site", values_to="LFC")%>%
  mutate(name=paste0("Phylum: ", Phylum, "", ", Genus:",Genus))
#plot this filtered data
infil<-ggplot(f_loc_com,aes(x=LFC, y=name, fill=site))+
  theme_minimal()+
  geom_col(position="dodge")+
  scale_fill_manual(values= palette)+
  labs(x="2-log fold change",y="")+
  geom_vline(xintercept = 0, size=2)+
  theme(text = element_text(size=14))
infil
ggsave(here("Figures/infil_deseq.png"),width = 7, height =7, units = "in")

treatment_ps<-subset_samples(genus, Timing == "After" & Treatment %in% c("NS", "WC"))
trtmnt_ds<-phyloseq_to_deseq2(treatment_ps,~Location + Treatment + Location:Treatment)
trtmnt_ds<-DESeq(trtmnt_ds)
#making normalized RA table
normalized_trt<-counts(trtmnt_ds,normalized=TRUE)%>%
  as_tibble(rownames = "ASV")%>%
  mutate_at(.,vars(-ASV),.funs = ~./sum(.))

tres1<-results(trtmnt_ds, contrast = c("Treatment", "WC", "NS"))
s_tres1<- tres1[which(tres1$padj<alpha),]
s_tres1<- cbind(as(s_tres1,"data.frame"),as(tax_table(genus)[rownames(s_tres1),],"matrix"))
s_tres1<- s_tres1[order(s_tres1$log2FoldChange),]
#impact treatment on ktya
tres2<-results(trtmnt_ds, contrast = list(c("Treatment_WC_vs_NS","LocationKTYA.TreatmentWC")))
s_tres2<- tres2[which(tres2$padj<alpha),]
s_tres2<- cbind(as(s_tres2,"data.frame"),as(tax_table(genus)[rownames(s_tres2),],"matrix"))
s_tres2<- s_tres2[order(s_tres2$log2FoldChange),] 
#impact treatment on HSP
tres3<-results(trtmnt_ds, contrast = list(c("Treatment_WC_vs_NS","LocationHSP.TreatmentWC")))
s_tres3<- tres3[which(tres3$padj<alpha),]
s_tres3<- cbind(as(s_tres3,"data.frame"),as(tax_table(genus)[rownames(s_tres3),],"matrix"))
s_tres3<- s_tres3[order(s_tres3$log2FoldChange),] 
  

#comparison from affect from 3 locations
loc_comp_trt <- tibble(ASV=row.names(s_tres1),LFC= s_tres1$log2FoldChange,Phylum=s_tres1$Phylum, Genus=s_tres1$Genus, loc="KTR")%>%
  rbind(tibble(ASV=row.names(s_tres2),LFC= s_tres2$log2FoldChange,Phylum=s_tres2$Phylum, Genus=s_tres2$Genus, loc="KTYA"))%>%
  rbind(tibble(ASV=row.names(s_tres3),LFC= s_tres3$log2FoldChange,Phylum=s_tres3$Phylum, Genus=s_tres3$Genus, loc="HSP"))
  
ggplot(loc_comp_trt, aes(x=LFC, y=Genus, fill=Phylum))+
  geom_col()+
  facet_grid(cols= vars(loc))+
  geom_vline(xintercept = 0,aes(size=2))

f_loc_com_trt<- pivot_wider(data = loc_comp_trt,id_cols = c("ASV", "Phylum","Genus"), names_from = loc, values_from = LFC)%>%
  subset(!is.na(KTR) & !is.na(HSP) &!is.na(KTYA))%>%
  pivot_longer(cols = c(KTYA,KTR,HSP), names_to= "site", values_to="LFC")%>%
  mutate(name=paste0("Phylum: ", Phylum, "", ", Genus:",Genus))
#try to look at RA of sign differential abudance ASVs 
#load metadata to be used to create avg
metadata<-readRDS(here("Objects/metadata.rds"))%>%
  as_tibble(rownames = "Sample")
trt_ASVs<-subset(f_loc_com_trt,duplicated(ASV)==FALSE)%>%
  select(ASV, Phylum, Genus)%>%
  inner_join(., normalized_trt, by="ASV")%>%
  #pivot so each sample has a value
  pivot_longer(cols = 4:57,names_to="Sample", values_to="rel_abund" )%>%
  #join to get metadata 
  inner_join(., metadata, by="Sample")


#plot this filtered data
ggplot(f_loc_com_trt,aes(x=LFC, y=name, fill=site))+
  geom_col(position="dodge")+
  theme_minimal()+
  scale_fill_manual(values= palette)+
  labs(x="2-log fold change",y="")+
  geom_vline(xintercept = 0, size=2)+
  theme(text = element_text(size=14))
ggsave(here("Figures/trt_deseq.png"), width = 7.5,height = 7.5, units = "in")


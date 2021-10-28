library(pacman)

p_load(tidyverse, ggplot2, here, readxl, ggpubr)
#renamed files to not have spaces *opened and saved(16S_plot30_0301202_Quantification_Cq_Results.xlsx)
#removed standard and "pos ctrl"
#NEED TO OPEN UP FILES AND SAVE
ref_16S <- read_excel(path= here("qPCR/16s_plot30_03162021_Quantification_Cq_Results.xlsx"),sheet = "0")%>%
  subset(!(Content %in% c("Std","Pos Ctrl","Neg Ctrl")))%>%
  select(1,5,7:9)%>%
  rename("Cq_16S"="Cq", "Cq_Mean_16S"="Cq Mean","std_dev_16S"="Cq Std. Dev")
nosZ <- read_excel(path = here("qPCR/nosZ_plot030_03162021_Quantification_Cq_Results.xlsx"),sheet = "0")%>%
  subset(!(Content %in% c("Std","Pos Ctrl","Neg Ctrl")))%>%
  select(1,5,7:9)%>%
  rename("Cq_nosZ"="Cq", "Cq_Mean_nosZ"="Cq Mean", "std_dev_nosZ"="Cq Std. Dev")

E_per_nosZ<- 57.4
E_per_16S<- 48.4
E_nosZ<- (E_per_nosZ/100)+1
E_16S<-(E_per_16S/100)+1
write_csv(qPCR, here("qPCR/qPCR_table.csv"))
qPCR<-inner_join(nosZ,ref_16S, by="Well")
saveRDS(qPCR,here("Objects/qPCR_30_raw.rds"))
#getting only one row per samples
sliceqPCR<- dplyr::slice(qPCR,1,4,7,10,13,16,19,22)
#normalizing to genomic DNA input, values from Qubit
dna_loads<-tibble(Sample=sliceqPCR$Sample.x, ng_DNA= c(8.55,7.95,9.8,7.3,3.06,2.49,5.05,4.56))
#saveRDS(dna_loads, here("qPCR/dna_loads"))
pfaffl.t<- read_csv(here("qPCR/qPCR_table.csv"))%>%
  mutate(., Timing=ifelse(Sample.x %in% c("kybn","kybw","ktbn","ktbw"), "before","after"))%>%
  mutate(Location=ifelse(Sample.x %in% c("kybn","kybw","kyan","kyaw"), "Kitayama","Kelly Thompson"))%>%
  mutate(Treatment=ifelse(Sample.x %in% c("kyaw","kybw","ktbw","ktaw"), "Woodchips","Native Soil"))
pfaffl<-read_csv(here("qPCR/qPCR_table.csv"))%>%
  select(., Sample.x, avg_ratio, sd_ratio)%>%
  subset(., !is.na(avg_ratio))
pfaffl<- mutate(pfaffl, Timing=ifelse(Sample.x %in% c("kybn","kybw","ktbn","ktbw"), "before","after"))%>%
  mutate(Location=ifelse(Sample.x %in% c("kybn","kybw","kyan","kyaw"), "Kitayama","Kelly Thompson"))%>%
  mutate(Treatment=ifelse(Sample.x %in% c("kyaw","kybw","ktbw","ktaw"), "Woodchips","Native Soil"))


t.test(ratio~Treatment, data=subset(pfaffl.t, Location =="Kitayama"& Timing=="after"))
t.test(ratio~Treatment, data=subset(pfaffl.t, Location =="Kelly Thompson"& Timing=="after"))
y<-expression(paste(italic("nosZ")," Ratio of After Infiltration to Before "))
#neither t-test was significant
palette<- readRDS(here("Objects/palette.RDS"))
ggplot(subset(pfaffl, Timing == "after"), aes(x=Treatment, y= avg_ratio, fill=Treatment))+
  geom_col()+
  facet_grid(cols = vars(Location))+
  theme_pubr(legend = "right")+
  geom_errorbar(aes(x=Treatment, ymin= avg_ratio-sd_ratio,ymax= avg_ratio+sd_ratio), width=0.2)+
  labs(x="", y= y)+
  scale_fill_manual(values = c(palette[4],palette[2]))

ggsave(here("Figures/qPCR_30cm_pfaffl.png"))

#changing nosZ Cq with no detection from zero to 40 (neg control) and adding informational columns about the samples
updated_qPCR <-sliceqPCR%>%
  #mutate(nosZCq_Mean=if_else(nosZCq_Mean==0, 40, nosZCq_Mean))%>%
  mutate(Sample.y=ifelse(Sample.y %in% c("kybn","kybw","ktbn","ktbw"), "before","after"))%>%
  mutate(Location=ifelse(Sample.x %in% c("kybn","kybw","kyan","kyaw"), "Kitayama","Kelly Thompson"))%>%
  mutate(Treatment=ifelse(Sample.x %in% c("kyaw","kybw","ktbw","ktaw"), "Woodchips","Native Soil"))%>%
  rename( "Timing"="Sample.y") %>%
  inner_join(., dna_loads, by= c("Sample.x"="Sample"))#%>%
 #mutate(Cq_Mean_nosZ= Cq_Mean_nosZ/ng_DNA)%>%
  #mutate(Cq_Mean_16S= Cq_Mean_16S/ng_DNA)

#making one file that can be manipulated for all calcs
fold_change<- pivot_wider(updated_qPCR, id_cols=c(Location,Treatment), 
                          names_from = Timing, values_from = c(Cq_Mean_nosZ, Cq_Mean_16S))%>%
  #dcq_nosZ is control(before) - sample(after)
  mutate(.,dcq_nosZ = Cq_Mean_nosZ_before - Cq_Mean_nosZ_after,
         dcq_16S = Cq_Mean_16S_before - Cq_Mean_16S_after,
         #dcq_after is sample(nosZ)-control(16S)
         dcq_after = Cq_Mean_nosZ_after - Cq_Mean_16S_after,
         dcq_before = Cq_Mean_nosZ_before - Cq_Mean_16S_before)%>%
  #ratio is dcq target(nosZ) / dcq ref(16S)
  mutate(ratio_nosZ = E_nosZ^dcq_nosZ ,
         ratio_16S = E_16S ^ dcq_16S,
         #ddcq is dcq sample(after) - dcq control (beofre)
         ddcq = dcq_after - dcq_before)%>%
  mutate(pfaffl_fold_change = signif(ratio_nosZ/ratio_16S,3),
         ddcq_fold_change = signif(2 ^ (-ddcq),3))

easy_to_read<-select(fold_change, 1:2,14:15)%>%
  pivot_longer(.,cols=c(pfaffl_fold_change,ddcq_fold_change), names_to="method",values_to="fold_change")%>%
  mutate(method= ifelse(method=="ddcq_fold_change", "Delta Delta Cq", "Pfaffl"))

#load palette

lab= expression(~Delta~Delta~Cq~Method)

ggplot(subset(easy_to_read, method=="Delta Delta Cq"), aes(y=fold_change,x= Treatment, fill=Treatment))+
  geom_col(width = .5, position = position_dodge(0.5))+
  facet_grid(cols=vars(Location))+
  geom_text(aes(label= fold_change), position = position_dodge(.5),vjust = -.3)+
  theme_pubr(legend = "right")+
  labs(x="", y= "2-Fold change")+
  scale_fill_manual(values = c(palette[4],palette[2]))+
  theme(text=element_text(size = 14))
ggsave(here("Figures/qPCR_30cm.png"))





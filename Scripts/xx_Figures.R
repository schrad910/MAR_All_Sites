library("tidyverse")
library("ggplot2")
library("ggpubr")
important_species <- subset(sigtabASV, !grepl("uncultured", sigtabASV$Species))%>%
  rownames_to_column(var = "featureID") %>%
  select(featureID) 

data_join <- metadata%>%
  select(SampleID, Location, Treatment, Timing, Treated, depth_below_plot_cm)
data_join$Timing[data_join$Timing == 0] <- "before"
data_join$Timing[data_join$Timing == 1] <- "after"
relabunplot<- filter(normalized_relabun,  normalized_relabun$feature.ID==important_species$featureID) %>% 
  select (-1)%>%
  pivot_longer(cols = (2:129), names_to = "SampleID", values_to = "rel_abun")%>%
  inner_join(.,data_join, by="SampleID") 
relabunplot[relabunplot == "Magnetospirillum_sp._enrichment_culture_clone_Van25"] <- "Magnetospirillum_sp."

ggplot(relabunplot, aes(x = rel_abun, y= Species, fill= Treatment )) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("Figures/Species_plot.png", width = 5, height = 5)
After_species <- subset(relabunplot, Timing =="after")
After_species <- subset(After_species, rel_abun < 0.03)

ggplot(After_species, aes(x = Species, y= rel_abun, fill= Treatment )) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Phylum plot

phylumplot<- pivot_longer(filter_Phycounts, cols = (2:129), names_to = "SampleID", values_to = "rel_abun")%>%
  inner_join(.,data_join, by="SampleID")
  saveRDS(phylumplot, file = "Objects/phylumplot.rds")

ggplot(phylumplot, aes(x = rel_abun, y= Phylum , fill= Timing )) +
  geom_boxplot() +
  theme_bw()+
  labs(title= "Phylum Relative Abundance Changes After Stimulated MAR",
       x="Relative Abundance")+
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        plot.caption = element_text(hjust = 0, size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))+
  scale_fill_manual(values= c("skyblue","skyblue4"))
ggsave("Figures/Phylum_plot.png", width = 5, height = 5)




After_phylum <- subset(phylumplot, Timing =="After")
ggplot(After_phylum, aes(x =rel_abun , y= Phylum , fill= Treatment )) +
  geom_boxplot()
ggsave("Figures/After_Phylum_plot.png", width = 5, height = 5)
 Proteo <- subset(After_phylum, Phylum == "Proteobacteria")
ggplot(Proteo, aes(x =rel_abun , y= Phylum , fill= Treatment )) +
   geom_boxplot()

##class plot
classplot_proteo <- pivot_longer(Classcounts, cols = (3:130), names_to = "SampleID", values_to = "rel_abun")%>%
  inner_join(.,data_join, by="SampleID")%>%
  subset(., Phylum == "Proteobacteria")
ggplot(classplot_proteo, aes(x =rel_abun , y= Class , fill= Treatment )) +
  geom_boxplot()
after_classplot_proteo<- subset(classplot_proteo, Timing == "after")
ggplot(after_classplot_proteo, aes(x =rel_abun , y= Class , fill= Treatment )) +
  geom_boxplot()
ggsave("Figures/After_Class_proteos.png", width = 5, height = 5)

classplot<- pivot_longer(Classcounts, cols = (3:130), names_to = "SampleID", values_to = "rel_abun")%>%
  inner_join(.,data_join, by="SampleID")

after_classplot<- subset(classplot, Timing == "after")
ggplot(after_classplot, aes(x =rel_abun , y= Class , fill= Treated)) +
  geom_boxplot()
ggsave("Figures/After_Class_treated.png", width = 5, height = 5)

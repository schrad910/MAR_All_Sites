library(tidyverse)
library(vegan)
library(ggvegan)
met_env<-read_tsv("~/Documents/SaltikovLab/Sequencing/All_Runs/metadata.txt")%>%
  column_to_rownames(.,var="SampleID")%>%
  select(Location,Treatment,Timing,C,N,Silt,Sand,Clay)
met_env$Treatment<- factor(met_env$Treatment, levels = c("NS","MX","WC","BC"), ordered = FALSE)
met_env$Timing <- factor(met_env$Timing, levels = c(0,1))

sigtabgen<- rownames_to_column(sigtabgen, var = "FeatureID")



linked_RA<- #inner_join(sigtabgen,normalized_relabun, by = c("FeatureID"="featureID"))%>%
  select(norm_tax_rel_abun, 13:141)%>%
  column_to_rownames(var="Genus")%>%
  t()

  cca<- cca( linked_RA~ C + N + Silt + Sand + Clay, data=met_env)

anova(cca)
anova(cca, by="axis")
anova(cca, by="terms")
anova(cca, by="margin")


fcca<- fortify(cca)
size<-1.8
ggplot(fcca, aes(x = CCA1, y = CCA2)) +
  geom_point(data = subset(fcca, Score == "species"),
             aes(shape = "species"), size = size) +
  geom_point(data = cbind(subset(fcca, Score == "sites"), Use = met_env$Treatment),
             aes(shape = "sites", color= Use), size = size) +
  scale_colour_brewer("Score", palette = "Set1") +
  coord_fixed() +
  theme(legend.position = "top")+
  geom_text(data = subset(fcca, Score == "biplot"), 
            aes(x=CCA1,y=CCA2,label=Label), size=2) + xlab("CCA1") + ylab("CCA2") +
  geom_segment(data = subset(fcca, Score == "biplot"),
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(1/2, 'picas')), colour = "red")
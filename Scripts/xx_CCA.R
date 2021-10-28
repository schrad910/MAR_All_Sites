library(tidyverse)
library(vegan)
library(ggvegan)
met_env<-read_tsv("~/Documents/SaltikovLab/Sequencing/MAR_All_Sites/Environment/metadata.txt")%>%
  #column_to_rownames(.,var="SampleID")%>%
  select(SampleID, Location,Treatment,Timing,C,N,Silt,Sand,Clay)

met_env$Treatment<- factor(met_env$Treatment, levels = c("NS","MX","WC","BC"), ordered = FALSE)
met_env$Timing <- factor(met_env$Timing, levels = c(0,1))





linked_RA<- #inner_join(sigtabgen,normalized_relabun, by = c("FeatureID"="featureID"))%>%
  select(norm_tax_rel_abun, 13:141)%>%
  column_to_rownames(var="Genus")%>%
  t()
filter_Classcounts<- subset(filter_Classcounts, Class != "Incertae_Sedis" | Class != 'uncultured')%>%
  select(-1)
  rownames(filter_Classcounts)<-NULL
t_classcount <-column_to_rownames(filter_Classcounts,var="Class")%>%
  t()
cca<- cca( t_classcount~ C + Sand + Clay, data=met_env)
KTR_phycounts<-t_phycounts[1:36,]
KTRenv<- subset(met_env, Location =="KTR")
KTR_cca<- cca( KTR_phycounts~ C + Sand + Clay, data=met_env)
anova(cca)
anova(cca, by="axis")
anova(cca, by="terms")
anova(cca, by="margin")
plot(cca, scaling = 3)
fcca<- fortify(cca)
size<-1.8
ggplot(fcca, aes(x = CCA1, y = CCA2)) +
  geom_point(data = subset(fcca, Score == "species"),
             aes(shape = "species"), size = size) +
  geom_point(data = cbind(subset(fcca, Score == "sites"), Use = met_env$Treatment),
             aes(shape = met_env$Location, color= Use), size = size) +
  scale_colour_manual(values=cal_palette("bigsur")) +
  coord_fixed() +
  theme(legend.position = "top")+
  geom_text(data = subset(fcca, Score == "biplot"), 
            aes(x=CCA1,y=CCA2,label=Label), size=2) + xlab("CCA1") + ylab("CCA2") +
  geom_segment(data = subset(fcca, Score == "biplot"),
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(1/2, 'picas')), colour = "red")

library(tidyverse)
library(vegan)
library(ggvegan)
library(here)
library(ggplot2)
library(ggpubr)



water_chem<- read_tsv(here("Environment/waterchem.txt"), 
                      col_type = cols(Timing = col_factor(levels = c("early","late")), 
                                  Plot = col_factor(levels = c("NS","MX","WC","BC")),
                                  Location = col_factor(levels = c("KTYA","KTR","HSP"))))%>%
                      rename(.,"TDS"="Total-Diss-Solids", "Alkalinity"="Alkalinity-CaCO3")


plot_water_chem<- pivot_longer(water_chem, cols= 6:25, names_to= "chem", values_to = "conc" )
                  
ggplot(plot_water_chem, aes(x=chem, y=conc, fill=Timing,))+
  geom_boxplot()+
  facet_grid(cols=vars(Location),
             rows = vars(Plot))+
  coord_flip()

Mn<-ggplot(water_chem, aes(x=Depth, y=Manganese, fill=Timing,))+
  geom_col(position = "dodge")+
  facet_grid(cols=vars(Location),
               rows = vars(Plot))

DOC<-ggplot(water_chem, aes(x=Depth, y=DOC, fill=Timing,))+
  geom_col(position = "dodge")+
  facet_grid(cols=vars(Location),
             rows = vars(Plot))
DOC
K<-ggplot(water_chem, aes(x=Depth, y=Potassium, fill=Timing,))+
  geom_col(position = "dodge")+
  facet_grid(cols=vars(Location),
             rows = vars(Plot))

EC<-ggplot(water_chem, aes(x=Depth, y=EC, fill=Timing,))+
  geom_col(position = "dodge")+
  facet_grid(cols=vars(Location),
             rows = vars(Plot))
Cr<-ggplot(water_chem, aes(x=Depth, y=Chromium, fill=Timing,))+
  geom_col(position = "dodge")+
  facet_grid(cols=vars(Location),
             rows = vars(Plot))
TDS<-ggplot(water_chem, aes(x=Depth, y=TDS, fill=Timing,))+
  geom_col(position = "dodge")+
  facet_grid(cols=vars(Location),
             rows = vars(Plot))
As<-ggplot(water_chem, aes(x=Depth, y=Arsenic, fill=Timing,))+
  geom_col(position = "dodge")+
  facet_grid(cols=vars(Location),
             rows = vars(Plot))
SO4<-ggplot(water_chem, aes(x=Depth, y=Sulfate, fill=Timing,))+
  geom_col(position = "dodge")+
  facet_grid(cols=vars(Location),
             rows = vars(Plot))
Fe<-ggplot(water_chem, aes(x=Depth, y=Iron, fill=Timing,))+
  geom_col(position = "dodge")+
  facet_grid(cols=vars(Location),
             rows = vars(Plot))
#incase Cr and DOC are important, taking out HSP incase its important
water_chem_noHSP<- subset(water_chem, Location !="HSP")%>%
  column_to_rownames(var = "Sample")
rda<- rda(water_chem_noHSP[,6:24])


water_chem_rda<- select(water_chem,1,6,7,9:23,25)%>%
  column_to_rownames(var = "Sample")
rda2<-rda(water_chem_rda)
palette<- readRDS(here("Objects/palette.RDS"))

#building a biplot! 
f_rda2<- fortify(rda2, axes=1:2)%>%
  left_join(.,water_chem[,1:5], by= c("Label"="Sample"))

take <- c('PC1', 'PC2') 
plot(rda2)
#Looked at the plot and made a list of the important species 
important_arrows<- c("EC", "TDS","Alkalinity", "Sulfate", "Magnesium","Calcium", "Manganese","Iron")
arrows <- f_rda2%>%
  subset( Score == 'species'& Label %in% important_arrows ) 
mul <- ggvegan:::arrowMul(arrows[, take],
                          subset(f_rda2, select = take, Score == 'sites'))
arrows[, take] <- arrows[, take] * mul

##plot showing the main drivers between KTYA& KTR
ggplot() +
  geom_point(data = subset(f_rda2, Score=="sites"),
             mapping = aes(x = PC1, y = PC2, col=Location, shape=Location)) + 
  geom_segment(data = arrows,
               mapping = aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text(data = arrows, # crudely push labels away arrow heads
            mapping = aes(label = Label, x = PC1 * 1.3, y = PC2 * 1.4)) +
  coord_fixed()+
  scale_color_manual(values = palette)+ 
  geom_abline(intercept = 0,slope = 0,linetype="dashed", size=0.8)+ 
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.8)+ 
  labs(x="PC1 (99.7%)", y="PC2 (0.3%)")+
  theme_pubr(legend = "right")
ggsave(here("Figures/waterchem_PCA.png"))

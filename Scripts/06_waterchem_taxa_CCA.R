## This script aims to connect water chemistry to the promotion of different species
## We only have water data for some of the after samples. Going to split into redox and physcial properties

library(pacman)
set.seed(3618)
p_load(vegan, tidyverse, ggvegan, ggplot2 ,ggpubr,here, phyloseq, ggrepel)
#load genus phyloseq
phylo<-readRDS(here("Objects/genera.rds"))

#looking at top 500 genus
top20gen = sort(tapply(taxa_sums(phylo), tax_table(phylo)[, "Genus"], sum), TRUE)[1:20]
gen20 = subset_taxa(phylo, Genus %in% names(top20gen))

species<- as(otu_table(gen20), "matrix")
tax<- as(tax_table(gen20), "matrix")%>%
  as_tibble(rownames = "Label")

#load env data
env<- read_tsv(here("Environment/metadata.txt"), 
         col_type = cols(Timing = col_factor(levels = c("Before","After")), 
                         Treatment = col_factor(levels = c("NS","MX","WC","BC")),
                         Location = col_factor(levels = c("HSP","KTR","KTYA"))))%>%
  select(-Sample1)
#subset to remove
env_waterchem<- subset(env,!is.na(DOC) & Treatment %in% c("NS", "WC"))%>%
  column_to_rownames(.,var="SampleID")
#env_waterchem2<-subset(env_waterchem, Location== "KTYA")%>%
  #select(.,-1:-13 )
  

species_waterchem<- subset(species, rownames(species) %in% rownames(env_waterchem))

## don't like  any of these models but may use later
#upr <- cca(species_waterchem ~ .,  data = env_waterchem2)
#lwr <- cca(species_waterchem ~ 1, data = env_waterchem2)

#mods <- ordistep(lwr, scope = formula(upr), trace = 0)
#stepwise selection give formula species_waterchem ~ Dissolved_solids + Na + CaCO3 + Zn + NO2 + As + Ca + pH + K + Cd
#72.5% of variance explained by constrained terms
#anova.cca(mods, by= "terms")
#vif.cca(mods)
#Diss_solids, NA, CaCO3, As, Ca, pH,K are all colinear.  removing #CaCO3
#mod1<- cca(species_waterchem~ Dissolved_solids+ Na + Zn +NO2 + As + Ca +pH + K + Cd, data= env_waterchem)
#anova.cca(mod1, by="terms") #Cd , K and As no longer sign
#mod2<-cca(species_waterchem~ Dissolved_solids+ Na + Zn +NO2  + Ca +pH , data= env_waterchem)
#anova.cca(mod1, by="terms") 
#vif.cca(mod2)
#autoplot(mod2)
#mod3<- cca(species_waterchem~ Dissolved_solids + Zn +NO2 +pH + Condition(Location) , data= env_waterchem)

#Makign model with redox- driven ions # including SO4 was codependent on NO3 
#model explains 50% 
redox_cca<- cca(species_waterchem~DOC + NO3 +  Mn + Fe , data = env_waterchem)
#anova.cca(redox_cca, by="terms")
#anova.cca(redox_cca, by="axis")
#anova.cca(redox_cca, by="margin")
vif.cca(redox_cca)#not codependent
summary(eigenvals(redox_cca))

f_redox_cca <- fortify(redox_cca, axes = 1:2, scaling= "species")%>%
  left_join(.,env[,1:7], by= c("Label"= "SampleID"))%>%
  left_join(.,tax, by="Label")
f_redox_cca1<-mutate(f_redox_cca, Genus= if_else(grepl('[0-9]', Genus), paste0(Family," ", Genus), Genus))

#make arrows
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(f_redox_cca, Score == 'biplot')  # take only biplot arrow scores
## multiplier for arrows to scale them to the plot range 
mul <- ggvegan:::arrowMul(arrows[, take],
                          subset(f_redox_cca, select = take, Score == 'species'))
arrows[, take] <- arrows[, take] * mul

#anova_test<-adonis2(species_waterchem~Treatment, data = env_waterchem, 
                 #   permutations = 10000 )

#anova_test

#plot with Phylums defined
palette<- colorRampPalette(readRDS(here("Objects/palette.RDS")))(9)
  palette1<- c("#89C0DB" , "#6d93a4" ,  "#5e3023" , "#3e5524" ,   "#edd9af" ,
             "#b3c5ab" , "#aa8b6f", "#799227" )


ggplot()+
  geom_text_repel(data = subset(f_redox_cca1, Score=="species" ),max.overlaps = 14,
         mapping = aes(x = CCA1, y = CCA2, label=Genus, col=Phylum), show.legend = F,fontface="italic", size=5) + 
  geom_point(data = subset(f_redox_cca, Score=="species" ),
             mapping = aes(x = CCA1, y = CCA2, color= Phylum))+
 #stat_ellipse(data=subset(f_redox_cca, Score=="sites"),geom="polygon", type="norm", alpha= 0.2,
  #                mapping = aes(x = CCA1, y = CCA2, col=Treatment, fill=Treatment),
       #       level=.8)+
  geom_segment(data = arrows,
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text(data = arrows, # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1*1.07, y = CCA2*1.05)) +
  #geom_point(data = subset(f_redox_cca, Score=="sites" & CCA2<5),
         # mapping= aes(x=CCA1, y=CCA2, shape=Location))+
  coord_fixed(ratio = 1.75)+
  scale_color_manual(values = palette1)+ 
  scale_fill_manual(values = palette1)+
  geom_abline(intercept = 0,slope = 0,linetype="dashed", size=0.8)+ 
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.8)+ 
  labs(x="CCA1(43.9%)", y="CCA2(3%)", caption = "ANOVA(model)<0.001")+
  theme_pubr(legend = "right")+
  theme(legend.text = element_text(face = "italic"), text=element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size=5)))


pstat <- permustats(anova(redox_cca))
summary(pstat)
ggsave(here("Figures/xx_redox_species.svg"), height = 10,width = 10, units = "in")







properties_cca<- cca(species_waterchem~K+ Cl+ pH  + Na , data = env_waterchem)
#anova.cca(properties_cca, by="terms")
#anova.cca(properties_cca, by="axis")
#anova.cca(properties_cca, by="margin")
#anova.cca(properties_cca)
vif.cca(properties_cca)
summary(eigenvals(properties_cca))
f_prop_cca <- fortify(properties_cca, axes = 1:2,)%>%
  left_join(.,env[,1:7], by= c("Label"= "SampleID"))%>%
  left_join(.,tax, by="Label")

f_prop_cca1<-mutate(f_prop_cca, Genus= if_else(grepl('[0-9]', Genus),Family, Genus))

#make arrows

p_arrows <- subset(f_prop_cca1, Score == 'biplot') %>%
  select(1:4)# take only biplot arrow scores
## multiplier for arrows to scale them to the plot range
p_mul <- ggvegan:::arrowMul(p_arrows[, take],
                          subset(f_prop_cca, select = take, Score == 'species'))
p_arrows[, take] <- p_arrows[, take] * p_mul






b<-ggplot()+
 geom_point(data = subset(f_prop_cca1, Score=="species"),
           mapping = aes(x = CCA1, y = CCA2, col=Phylum)) + 
  geom_text(data = subset(f_prop_cca1, Score=="species" ),
                mapping = aes(x = CCA1, y = CCA2, label=Genus, col=Phylum), show.legend = F) +
  #stat_ellipse(data=subset(f_prop_cca, Score=="species"),
               #mapping = aes(x = CCA1, y = CCA2, col=Phylum))+
  geom_segment(data = p_arrows,
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text(data = p_arrows, # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1*1.1, y = CCA2*1.1)) +
  #geom_point(data = subset(f_prop_cca, Score=="sites" & CCA2<5),
   #          mapping= aes(x=CCA1, y=CCA2, col= Treatment, shape=Location))+
  #stat_ellipse(data=subset(f_prop_cca, Score=="sites"),
  #mapping = aes(x = CCA1, y = CCA2, col=Treatment),type = "t")+
  coord_fixed()+
  scale_color_manual(values = palette)+ 
  geom_abline(intercept = 0,slope = 0,linetype="dashed", size=0.8)+ 
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.8)+ 
  labs(x="CCA1 (44%)", y="CCA2 (7.8%)")+
  theme_pubr(legend = "right")

ggsave(here("Figures/xx_CCA_waterchem_sup.png"))





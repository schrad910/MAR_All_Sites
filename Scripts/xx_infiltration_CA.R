### This script is trying to view how infiltration drives bacterial communities.... I need to 
### figure out how exactly to test this given that i only have C, N, and GS for before/after samples
### I think that i shouldn't be doing CCA and maybe do NDMS or something like CA
### Also need to look into doing permutations tests with plots to separate locations. 
library(pacman)
p_load(tidyverse, vegan, ggvegan, here, ggplot2, ggpubr, phyloseq)

set.seed(3618)
#loading the glommed phylum phyloseq object
phylo<- readRDS(here("Objects/genera.rds"))

#pulling the table from phylo cause it combines the genuses and looks nicer
table<- as(otu_table(phylo), "matrix")
tax<- as(tax_table(phylo), "matrix")%>%
  as_tibble(rownames = "Label")
#brining enviro data adn subsetting to only include NS and WC
env<- read_tsv(here("Environment/metadata.txt")) %>%
  select(-Sample1)


         

#subsetting to only include NS samples(aka look at infil)
env_ns<- subset(env, Treatment =="NS")
table_ns<-subset(table, rownames(table) %in% env_ns$SampleID)%>%
  as_tibble(rownames = NA)
#adding in soil differences to constrain plot
ca<- cca(table_ns, data=env_ns)

 # fortify the ordinations
f_ca<-fortify(ca, axes = 1:2, scaling=2)%>%
  left_join(.,tax, by="Label")%>%
  left_join(., env_ns, by= c("Label"= "SampleID"))
phylum<-tax_glom(phylo, taxrank = "Phylum")%>%
  subset_taxa(., Phylum=="Proteobacteria" | Phylum =="Actinobacteriota")

phylum1<- as(otu_table(phylum), "matrix")%>%
  as_tibble(rownames = "SampleID")%>%
  rename(.,"Actinobacteriota" ="ASV22", "Proteobacteria"="ASV41")%>%
  subset(., SampleID %in% env_ns$SampleID)%>%
  left_join(.,env_ns, by  = "SampleID")%>%
  column_to_rownames(., var="SampleID")
  
community_shift<- envfit(ca~Proteobacteria + Actinobacteriota, data=phylum1, choices=1:2, scaling=1, permutations=1000)
arrow<-as.data.frame(scores(community_shift, "vectors"))* ordiArrowMul(community_shift)


#anova test tells us that timing isn't significantly differnt
anova_test<-adonis2(table_ns~Timing, data = phylum1, permutations = 10000 )
anova_test
#Before and after are not significantly different
palette<- readRDS(here("Objects/palette.RDS"))

ggplot() +
  geom_point(data = subset(f_ca, Score=="sites"),
             mapping = aes(x = CA1, y =CA2,  shape=Timing, col=Timing), size=2) + 
  stat_ellipse(data= subset(f_ca, Score =="sites"), mapping = aes(x = CA1, y = CA2, col=Timing, fill=Timing),
            geom="polygon", type="norm", alpha= 0.2, show.legend = F, level = .7)+
  #stat_ellipse(data= subset(f_ca, Phylum =="Proteobacteria" | Phylum =="Actinobacteriota"), mapping = aes(x = CA1, y = CA2, col=Phylum, fill=Phylum),
   #            geom="polygon", type="norm", alpha= 0.2, show.legend = F, level = .8)+
  geom_point(data = subset(f_ca, Score=="species"),
            mapping= aes(x=CA1, y=CA2), shape=3, size=.2)+
  coord_fixed(ratio= .8)+
  scale_color_manual(values = c("#593F28", "#89C0DB"))+ 
  scale_fill_manual(values = c("#593F28", "#89C0DB"))+ 
  geom_abline(intercept = 0,slope = 0,linetype="dashed", size=0.8)+ 
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.8)+ 
  labs(x="CA1 (21.7%)", y="CA2 (12.7%)", caption = "ANOVA(Timing)< 0.1" )+
  theme_pubr(legend = "right")+
  geom_text(data = arrow, aes(x = CA1*.55 , y = CA2*.55), 
            label = row.names(arrow), colour = "black", 
  ) + 
  geom_segment(aes(x = 0, y = 0, xend = CA1*.5, yend = CA2*.5), 
               data = arrow, alpha=.7,  colour = "black",
               arrow = arrow(length = unit(0.01, "npc")))

##NMDS is not more informative at all.... 
nmds<- metaMDS(table_ns)

f_nmds<-fortify(nmds)
f_sites_nmds<- subset(f_nmds, Score=="sites")
f_sites_nmds$Label<-rownames(table_ns)

f_nmds<-subset(f_nmds,Score=="species")%>%
 rbind(f_sites_nmds)%>%
  left_join(.,env_ns, by= c("Label"= "SampleID"))%>%
  left_join(.,tax[,1:3], by="Label")

ggplot() +
  geom_point(data = subset(f_nmds, Score=="sites"),
             mapping = aes(x = NMDS1, y =NMDS2,  shape=Timing, col=Location)) + 
  #stat_ellipse(data= subset(f_nmds, Score =="sites"), mapping = aes(x = NMDS1, y = NMDS2, col=Timing, fill=Timing),
   #            geom="polygon", type="norm", alpha= 0.2, show.legend = F, level = .8)
#stat_ellipse(data= subset(f_nmds, Phylum =="Proteobacteria" | Phylum =="Actinobacteriota"), mapping = aes(x = NMDS1, y = NMDS2, col=Phylum, fill=Phylum),
 #            geom="polygon", type="norm", alpha= 0.2, show.legend = F, level = .8)
  #geom_point(data = subset(f_nmds, Score=="species"),
   #          mapping= aes(x=NMDS1, y=NMDS2), shape=3, size=.2)+
  coord_fixed(ratio= 1.4)
  #scale_color_manual(values = palette)+ 
  #scale_fill_manual(values = palette)+ 
  #geom_abline(intercept = 0,slope = 0,linetype="dashed", size=0.8)+ 
  #geom_vline(aes(xintercept=0), linetype="dashed", size=0.8)+ 
  #labs(x="CA1 (*fill)", y="CA2 (*fill)")+
  #theme_pubr(legend = "right")



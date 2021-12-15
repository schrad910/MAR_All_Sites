library(pacman)

p_load(vegan, tidyverse, ggvegan, ggplot2,here, phyloseq, ggpubr, ggrepel)
set.seed(3618)

phylo<-readRDS(here("Objects/genera.rds"))

before<-subset_samples(phylo, Timing == "Before"& Treatment %in% c("NS", "WC"))



reads_before<-as(otu_table(before), "matrix")
taxa_before<-as(tax_table(before), "matrix")%>%
  as_tibble(rownames = "Label")
env<- read_tsv(here("Environment/metadata.txt"), 
               col_type = cols(Timing = col_factor(levels = c("Before","After")), 
                               Treatment = col_factor(levels = c("NS","MX","WC","BC")),
                               Location = col_factor(levels = c("HSP", "KTR","KTYA"))))%>%
  select(-Sample1)

env<-subset (env, Timing=="Before" & Treatment %in% c("NS", "WC")) 



test<-adonis2(reads_before~Location, data = env, permutations = 10000 )
#location is significant
#making model with all inputs sand and clat added first becuase that is waht stepwise chose
#the rest are in order of differences between averages. 
model<- cca(reads_before~Sand + Clay + Silt + C +N , data=env )
#model= reads_before ~Clay + Sand
#anova.cca(model) #significant model
#anova.cca(model, by="axis") #top2 axis are significant
#anova.cca(model, by="terms") #sand and clay <.001, c < .2
#vif.cca(model)##values all <10 so  no collinearity

summary(eigenvals(model)) #constrained explains 31.2% of variance
#RsquareAdj(model) #adjusted is 20.5%



f_model<-fortify(model, axes = 1:2, scaling= "sites")%>%
  left_join(.,env[,1:7], by= c("Label"= "SampleID"))%>%
  left_join(.,taxa_before, by="Label")
#make arrows
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(f_model, Score == 'biplot')  # take only biplot arrow scores
## multiplier for arrows to scale them to the plot range 
mul <- ggvegan:::arrowMul(arrows[, take],
                          subset(f_model, select = take, Score == 'sites'))
arrows[, take] <- arrows[, take] * mul
palette<- readRDS(here("Objects/palette.RDS"))

ggplot()+
  #geom_point(data = subset(f_model, Score=="species" ),
   #          mapping = aes(x = CCA1, y = CCA2), shape=3) + 
  stat_ellipse(data=subset(f_model, Score=="sites"),geom="polygon", type="t", alpha= 0.2,
              mapping = aes(x = CCA1, y = CCA2, col=Location, fill=Location, linetype=Location), level=.95)+
  geom_segment(data = arrows,
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text(data = arrows, # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1*1.15, y = CCA2*1.15) ) +
  geom_point(data = subset(f_model, Score=="sites" & CCA2<5),
             mapping= aes(x=CCA1, y=CCA2, col= Location, shape= Location))+
  coord_fixed(ratio=1.2)+
  scale_color_manual(values = palette)+ 
  scale_fill_manual(values = palette)+
  scale_shape_manual(values = c(15,17,19))+
  geom_abline(intercept = 0,slope = 0,linetype="dashed", size=0.8)+ 
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.8)+ 
  labs(x="CCA1 (12.8%)", y="CCA2 (11.2%)" , caption = "PERMANOVA(Location)<.001")+ 
  theme_pubr(legend = "right")

ggsave(here("Figures/02_before_cca.svg"))


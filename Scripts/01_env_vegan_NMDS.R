library(pacman)

p_load(vegan, tidyverse, ggvegan, ggplot2,here, ggpubr)
set.seed(3618)
env<-read_tsv(here("Environment/metadata.txt"))%>%
  as_tibble()%>%
  #reducing # of samples 
  subset(Replicate==1 & Timing== "After")%>%
  column_to_rownames(var = "SampleID")


env_PCA<-rda(env[,9:13], scale=T)
env_NMDS<-metaMDS(env[,9:13])
vectors<-envfit(env_NMDS~C+N+Sand+Silt+Clay, data=env, choices=1:2, scaling=1, permutations=1000)
arrow<-as.data.frame(scores(vectors, "vectors"))* ordiArrowMul(vectors)
palette<- readRDS(here("Objects/palette.RDS"))

#building a biplot! 
#env_data<- fortify(env_PCA, axes=1:2, display = "sites")
f_env<-fortify(env_PCA, axes=1:2)%>%
  inner_join(., rownames_to_column(env, var="Sample"), by= c("Label"="Sample"))

#test to make sure differences in location are significant.
anova_test<-adonis2(env[,9:13]~Location, data = env, permutations = 10000 )
anova_test

ggplot() +
  geom_point(data = f_env,
             mapping = aes(x = PC1, y = PC2, col=Location, shape=Location,
                           fill=Location)) + 
  coord_fixed(ratio=.5)+
  scale_color_manual(values = palette)+
  scale_fill_manual(values = palette)+
  geom_abline(intercept = 0,slope = 0,linetype=2, size=0.8, colour= 'grey')  + 
  geom_vline(aes(xintercept=0), linetype=2, size=0.8, colour= 'grey')+ 
  labs(x="NMDS1", y="NMDS2", caption="ANOVA(Location) < 0.001")+
  scale_shape_manual(values = c(15,17,19))+
  theme_pubr(legend = "right")+
  geom_text(data = arrow, aes(x = NMDS1*1.1 , y = NMDS2*1.1), 
            label = row.names(arrow), colour = "black", 
            ) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = arrow, alpha=.7,  colour = "black",
               arrow = arrow(length = unit(0.01, "npc")))+ 
  stat_ellipse(data= f_env, aes(x=PC1, y=PC2, col=Location, fill=Location, linetype=Location), 
               geom="polygon", type="t", alpha= 0.2, level = .95)



 ggsave(here("Figures/01_PCA_env.pdf"))
  
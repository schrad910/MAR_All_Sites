library(vegan)
library(ggvegan)
library(tidyverse)
library(ggplot2)
library(gridExtra)
theme_light()
#making only 1 replicate so it is less busy
metadata<-readRDS("Objects/metadata.rds")
metadata_env <- subset(metadata, metadata$Replicate == 1)
#changing depth to a factor 
breaks <- c(0,40,70,120)
names<- c("10-39 cm", "40-69 cm", "70-110 cm")
metadata_env$Depth<- cut(metadata_env$depth_below_plot_cm, breaks=breaks,labels = names)

env_1 <- ggplot(metadata_env, aes(x=C, y=N, color= Location, shape=Depth)) + 
  geom_jitter() + 
  scale_color_manual(values=c("#93AFBB","#E8CBAE","#21554F"))+
  theme_light()+
  labs(x="% weight C", y="% weight N") +
  theme(legend.position="none")
env_2 <- ggplot(metadata_env, aes(x=Sand, y= Silt + Clay, color= Location, shape=Depth)) + 
  geom_jitter() + 
   scale_color_manual(values=c("#93AFBB","#E8CBAE","#21554F"))+
  theme_light()+
  labs(x= "% Sand", y = "% Silt and Clay")+
  theme(legend.position="none")

a<-arrangeGrob(env_1,env_2, ncol=2)
#legend saved separately
ggsave("Figures/environmental_data.png", a,width = 8, height = 5)

#making ndms 
library(vegan)
env_data <- column_to_rownames(metadata_env, "SampleID")
env_data<- select(env_data, 10:16)
  env_data <- select(env_data, -Depth)
rda<- rda(env_data)                    
biplot(rda)                    

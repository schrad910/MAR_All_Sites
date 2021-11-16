library(pacman)
set.seed(3618)
p_load( tidyverse,  ggplot2,here, ggpubr)

env<- read_tsv(here("Environment/metadata.txt"))%>%
  subset(., Replicate == "1" )%>%
  select(Location, Treatment, Timing, depth_below_plot_cm, C, N, Clay, Silt, Sand)%>%
  subset(.,Treatment %in% c("NS","WC"))
env<-mutate(env, Treatment= ifelse(Treatment == "NS", "Native Soil plot", "Wood chip amended plot"))
avg<- env%>%
  group_by(Location)%>%
  mutate(C_avg= mean(C), N_avg= mean(N), sand_a= mean(Sand), silt_a=mean(Silt), clay_a=mean(Clay))
palette<- readRDS(here("Objects/palette.RDS"))
##making plots looking at changes between C and N due to iniliftratin and location.
C_plot<-ggplot()+
  geom_boxplot(data= env, aes(x=Location, y= C, col=Timing), position = "dodge")+
  coord_fixed(3)+
  theme_pubr(legend = "right")+
  labs(y= "C (%-weight)")+
  scale_color_manual(values=palette)
N_plot<-ggplot()+
  geom_boxplot(data= env, aes(x=Location, y= N, col=Timing), position = "dodge")+
  coord_fixed(33)+
  theme_pubr(legend = "right")+
  labs(y= "N (%-weight)")+
  scale_color_manual(values=palette)
ggarrange(C_plot,N_plot, labels = c("A","B"), ncol=2, nrow = 1)
ggsave(here("Figures/SupplementaryFigs/SupFig02.tiff"), height= 2, width = 7.5, units = "in")


ggplot()+
  geom_line(data = env, aes(x=depth_below_plot_cm,  y=Sand, col="Sand"), linetype=2)+
  geom_point(data = env, aes(x=depth_below_plot_cm,  y=Sand, col="Sand", shape="Sand"))+
  geom_line(data = env, aes(x=depth_below_plot_cm,  y=Silt, col="Silt"), linetype=1)+
  geom_point(data = env, aes(x=depth_below_plot_cm,  y=Silt, col="Silt", shape= "Silt"))+
  geom_line(data = env, aes(x=depth_below_plot_cm,  y=Clay, col="Clay"), linetype=4)+
  geom_point(data = env, aes(x=depth_below_plot_cm,  y=Clay, col="Clay", shape= "Clay"))+
  scale_x_reverse() +
  coord_flip()+ 
  theme_pubr(legend = "right")+
  facet_grid(cols = vars(Location), rows= vars(Treatment))+
  scale_color_manual(name= "Grain Size", labels= c("Clay", "Sand", "Silt"),values= c("#593F28", "#D0C1B0" ,"#739D6E"))+
  scale_shape_manual(name= "Grain Size", labels= c("Clay", "Sand", "Silt"), values = c(15,17,19))+
  labs(y="Percent", x="Depth below plot (cm)")
ggsave(here("Figures/SupplementaryFigs/SupFig03.tiff"))  

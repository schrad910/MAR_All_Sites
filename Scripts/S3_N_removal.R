library(pacman)
p_load(tidyverse, readxl, ggplot2, here, ggpubr, scales)

native_soil <- read_xlsx(path = here('Environment/NRemovalComparison.xlsx'), sheet = 'NS', range = "G3:L18", 
                 col_names = c("N_load_HSP","Frac_rem_HSP","N_load_KTR","Frac_rem_KTR","N_load_KTYA","Frac_rem_KTYA")) %>%
  mutate(., sample= paste0("Sample",seq(nrow(.))))%>%
  pivot_longer(cols = 1:6, names_to= c("Measurement", "Location"),
               names_pattern= "(.*_.*)_(.*)",
               values_drop_na=T, 
               values_to="value")%>%
  pivot_wider(., id_cols= c("Location","sample"), names_from="Measurement", values_from ="value")%>%
  select(-sample)
load_title<-expression("Initial N load (g-N/m"^2~d*")")

woodchips<- read_xlsx(path = here('Environment/NRemovalComparison.xlsx'), sheet = 'WC', range = "G3:L19", 
                 col_names = c("N_load_HSP","Frac_rem_HSP","N_load_KTR","Frac_rem_KTR","N_load_KTYA","Frac_rem_KTYA"))%>%
  mutate(., sample= paste0("Sample",seq(nrow(.))))%>%
  pivot_longer(cols = 1:6, names_to= c("Measurement", "Location"),
               names_pattern= "(.*_.*)_(.*)",
               values_drop_na=T, 
               values_to="value")%>%
  pivot_wider(., id_cols= c("Location","sample"), names_from="Measurement", values_from ="value")%>%
  select(-sample)

palette<-readRDS(here("Objects/palette.RDS"))

ggplot()+
  geom_point(data= native_soil, aes(x=N_load, y=Frac_rem, shape= Location, col="Native Soil")) +
  geom_smooth(data= native_soil,method="lm",se=FALSE, aes(x=N_load, y=Frac_rem,  col="Native Soil"))+
  #stat_regline_equation(data = native_soil, aes(x=N_load, y=Frac_rem,  col="Native Soil",
   #                     label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~")), label.y.npc = "bottom", show.legend = F)+
  #stat_cor(data = native_soil, aes(x=N_load, y=Frac_rem,  col="Native Soil"), 
#label.y.npc = "bottom")+
  geom_point(data= woodchips, aes(x=N_load, y=Frac_rem, shape= Location, col="Woodchips"))+
  geom_smooth(data= woodchips,method="lm", se=FALSE,  aes(x=N_load, y=Frac_rem, col="Woodchips"))+
  #stat_regline_equation(data= woodchips,aes(x=N_load, y=Frac_rem, col="Woodchips",
   #                                         label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~")), show.legend = F)+
  #stat_cor(data= woodchips,aes(x=N_load, y=Frac_rem, col="Woodchips"))+
  scale_color_manual(name = "Treatment", 
        values = c("Native Soil" = "#D0C1B0" , "Woodchips" = "#80654C"))+
  geom_abline(intercept = 0,slope = 0,linetype=2, size=0.8, colour= 'grey')  + 
  labs(x=load_title, y="Fraction N removed")+
  scale_shape_manual(values = c(15,17,19))+
  theme_pubr(legend = "right")+ 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), 
                limits= c(10^0.6,10^2.2))+
  annotation_logticks(sides="b")+ 
  coord_fixed(ratio=1.1)



  ggsave(here("Figures/N_rem.png"))
total<-rbind(mutate(native_soil, treatment= "Native soil"), mutate(woodchips, treatment="Woodchips"))

             
  
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
  select(-sample)%>%
  mutate(Rate= case_when(N_load<=10 ~ "Slow",
                         N_load<= (10^1.5) ~ "Medium",
                         N_load>(10^1.5) ~"Fast"))
load_title<-expression("Initial N load (g-N/m"^2~d*")")

woodchips<- read_xlsx(path = here('Environment/NRemovalComparison.xlsx'), sheet = 'WC', range = "G3:L19", 
                 col_names = c("N_load_HSP","Frac_rem_HSP","N_load_KTR","Frac_rem_KTR","N_load_KTYA","Frac_rem_KTYA"))%>%
  mutate(., sample= paste0("Sample",seq(nrow(.))))%>%
  pivot_longer(cols = 1:6, names_to= c("Measurement", "Location"),
               names_pattern= "(.*_.*)_(.*)",
               values_drop_na=T, 
               values_to="value")%>%
  pivot_wider(., id_cols= c("Location","sample"), names_from="Measurement", values_from ="value")%>%
  select(-sample)%>%
  mutate(Rate= case_when(N_load<=10 ~ "Slow",
                         N_load<= (10^1.5) ~ "Medium",
                         N_load>(10^1.5) ~"Fast"))

palette<-readRDS(here("Objects/palette.RDS"))

ggplot()+
  geom_point(data= native_soil, aes(x=N_load, y=Frac_rem, shape= Location, col="Native Soil")) +
  geom_smooth(data= native_soil,method="gam",se=FALSE, aes(x=N_load, y=Frac_rem,  col="Native Soil"))+
  geom_point(data= woodchips, aes(x=N_load, y=Frac_rem, shape= Location, col="Woodchips"))+
  geom_smooth(data= woodchips,method="gam", se=FALSE,  aes(x=N_load, y=Frac_rem, col="Woodchips"))+
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



#ggsave(here("Figures/S3_N_rem.pdf"))
total<-rbind(mutate(native_soil, Treatment= "Native Soil"), mutate(woodchips, Treatment="Woodchips"))
facet_labels <- function(variable, value) {
  labels <- as.character(value)
  labels[labels == '(all)'] <- 'All'
  return (labels)
}
ggplot(data=total,
       aes( x= Treatment, y=Frac_rem))+
  geom_boxplot(position="dodge",outlier.shape = NA)+ 
  geom_point(aes(color=Treatment, shape=Treatment), position=position_jitterdodge())+
  facet_grid(cols=vars(Location), margins=T, labeller = facet_labels)+
  geom_hline(aes(yintercept=0), linetype=2, color= "gray")+
  theme_pubr(legend = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = NA, color = "black"))+
  scale_y_continuous(breaks = c(-.6,-0.4, -0.2,0, .2,.4,.6) )+
  stat_compare_means(aes( x=Treatment, y=Frac_rem), 
                  method= "t.test",label="p.signif",label.x = 1.5, label.y = .6,
                   hide.ns = T)+
  scale_color_manual(values = c("#D0C1B0","#80654C"))+
  #scale_fill_manual(values = c("#D0C1B0","#80654C"))+
  labs(x="Treatment",y="Fraction of Dissolved Nitrogen Removed")+
  coord_fixed(ratio=3.8)

ggsave(here("Figures/Round2/Figure01.tiff"))

novo_counts<- filter(normalized_relabun, Genus=="Novosphingobium"| Genus=="Sphingobium")%>%
  select(-feature.ID, -Domain, -Phylum, -Family, -Class, -Order, -Species)%>%
  pivot_longer(!Genus,names_to="sample", values_to="rel_abun" )%>%
  aggregate(rel_abun~sample+Genus, data=., sum)%>%
  inner_join(., met_env, by=c("sample"="SampleID"))%>%
  filter(rel_abun > 0.005)


ggplot(novo_counts, aes(x=sample,y=rel_abun, fill= Genus))+
  geom_col(position = "dodge")+
  theme(axis.text.x = element_text(angle = 90))

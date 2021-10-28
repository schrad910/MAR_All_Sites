library("tidyverse")
library("ggplot2")
library("ggpubr")
#phylum order
#sorts descendingy by log2fold change

#genus order
x= tapply(sigtabASV$log2FoldChange, sigtabASV$Species, function(x) max(x))%>%
  sort(.,TRUE)
sigtabASV$Species = factor(as.character(sigtabASV$Species), levels = names(x))

theme_classic()
ggplot(sigtabASV, aes(x=log2FoldChange, y= Species, color= Phylum)) +
  geom_vline(xintercept=0, color= "black", size=.5)+
  geom_point(size=4)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5))

inner_join()
ggplot(, aes(x = Species, y= rel_abun , fill= Treatment )) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90,  vjust= 0.2, hjust=.4))

#look at just proteos
proteo = subset(sigtabASV, Phylum == "Proteobacteria")
x=tapply(proteo$log2FoldChange, proteo$Species, function(x) max(x))%>%
  sort(., TRUE)



ggplot(proteo, aes(x=log2FoldChange, y= Species, color= Class)) +
  geom_vline(xintercept=0, color= "black", size=.5)+
  geom_point(size=4)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5))


#prevelancedf = apply(X = otu_table(phylo),
 #                    MARGIN = 1,
  #                   FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
#prevelancedf = data.frame(Prevalence = prevelancedf,
 #                         TotalAbundance = taxa_sums(phylo),
  #                        tax_table(phylo))
library("tidyverse")
#phylum order
#sorts descendingy by log2fold change
sort(.,TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
#genus order
x= tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))%>%
  sort(.,TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))

ggplot(sigtabgen, aes(x=log2FoldChange, y=Genus, color= Phylum)) +
  geom_vline(xintercept=0, color= "black", size=.5)+
  geom_point(size=1)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5))

#look at just proteos
proteo = subset(sigtabgen, Phylum == "Proteobacteria")
x=tapply(proteo$log2FoldChange, proteo$Species, function(x) max(x))%>%
  sort(., TRUE)
proteo$Phylum = factor(as.character(proteo$Phylum), levels = names(x))



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
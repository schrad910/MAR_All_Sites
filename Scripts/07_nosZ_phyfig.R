library(pacman)
p_load(tidyverse,ggplot2, phyloseq, vegan, here, ggvegan, ggtree, tidytree, ape, ggtext)
set.seed(3618)


genus<- readRDS(here("Objects/nosZ_phylo.rds"))
#fit<- readRDS(here("Objects/nosz_phytree.RDS"))

rel_genus <- transform_sample_counts(genus, function(otu) otu/sum(otu))


#g_table<- as(otu_table(genus), "matrix")
#taxa<- as(tax_table(genus), "matrix")%>%
  #as.data.frame()

#metadata<- as(sample_data(genus),"data.frame")

newick<- read.tree(here("nosZ/newick.nw")) 
ref_taxa<-read_csv(here("nosZ/ref_taxa.csv"))
tree<- fortify(newick)
tree$label<-gsub("_", " ", tree$label)

  


table<-fortify(rel_genus)%>%
  subset(., isTip==TRUE)%>%
  select(-x,-y,-parent,-node,-isTip, -branch.length, -angle, -branch)
joined_tree<-right_join(table, tree, by="label")
#adding in taxa 
joined_tree[122:132,8:13]<-ref_taxa[,2:7]


joined_tree$Genus<-with(joined_tree, ifelse(is.na(Genus)& isTip==TRUE, "Unclassified", 
                                    paste0("*", Genus, "*")))
joined_tree$Family<-with(joined_tree, ifelse(is.na(Family)& isTip==TRUE, "Unclassified", 
                                            paste0("*", Family, "*")))


palette<-c("#ccc0b6","#564b40", "#355e3b", "#88a47d" , "#89C0DB",
            "#b38b6d", "#8c8c8c", "#80654C","#c2d1bc", "#000000" )



p<-ggtree(joined_tree, ladderize = F)+
  geom_tiplab( aes(label = label, color=Family), 
              align=T,size=5,show.legend=F , offset = .01)+
  scale_color_manual(values=palette)+
  geom_nodelab(aes(label=label), nudge_x = -.04, nudge_y = .42)+
  geom_tippoint(aes(color=Family), size=4, )+
  geom_treescale(x=0, y=-0.5, width=.1)+
  xlim_tree(1.5)+
  theme(legend.position = "left",legend.text = element_markdown(size = 8))+ 
  guides(color = guide_legend(override.aes = list(size=5)))
  heatmap<- mutate(joined_tree, Name= paste(Location, Timing, Treatment, sep = "; "))%>%
    subset(., isTip==TRUE)%>%
    select(., label, Name, Abundance) %>%
    subset(., !is.na(label))%>%
    pivot_wider(., id_cols= "label", names_from ="Name", values_from= "Abundance")%>%
    column_to_rownames(var="label")%>%
    select(-9)
  
  
  #change samples from NA to 0 for formatting purposes
  heatmap[1:17, 1:8][is.na(heatmap[1:17, 1:8])] <- 0
  #reorder columns
  heatmap<- select(heatmap, 'KTR; Before; NS', 'KTR; After; NS','KTR; Before; WC','KTR; After; WC', 
                   'KTYA; Before; NS', 'KTYA; After; NS', 'KTYA; Before; WC','KTYA; After; WC')

gheatmap(p,heatmap,
           offset=.52,
           width=0.3,
           color = "grey90",
           font.size = 4,
           legend_title="Abundance",
           colnames_angle=90,
           hjust=1) +
  scale_fill_gradient(low =  "#89C0DB" ,
                        high = "#593F28", na.value = "white", name="Relative Abundance")+
  ggplot2::ylim(-5, 30) +
  theme(legend.position = "right", legend.text = element_markdown(size=12) )
  
ggsave(here("Figures/07_nosZ_phylo_heatmap.svg"), width= 14, height = 10, units = "in")  


#library(ANCOMBC) ##neither treatment or timign is significant
#out = ancombc(phyloseq= genus, formula = "Timing", p_adj_method = "holm", zero_cut=.9,
#              struc_zero=FALSE,   alpha=0.05) 

#out_df <- data.frame(
#  Species = row.names(out$res$beta),
#  coef = unlist(out$res$beta),
 # se = unlist(out$res$se),
#  W = unlist(out$res$W),
#  p_val = unlist(out$res$p_val),
#  q_val = unlist(out$res$q_val),
#  diff_abn = unlist(out$res$diff_abn))

#phylo<- readRDS("Objects/genera.rds")
#phylo<-transform_sample_counts(phylo, function(otu) otu/sum(otu))

#table_16S<-as(otu_table(phylo), "matrix")%>%
#  t()
#taxa_16S<-as(tax_table(phylo), "matrix")%>%
#  as_tibble(rownames = "asv") 

#s_taxa<- subset(taxa_16S, Genus %in% joined_tree$Genus)
#asv<- s_taxa$asv

#table_16S_s<-subset(table_16S, rownames(table_16S) %in% asv)

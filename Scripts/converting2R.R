#This script pulls in files made in qiime and sets them to be used by R
library(qiime2R)
library(tidyverse)
library(phyloseq)

#renaming columns with problematic names
metadata<-read_tsv("metadata.txt", col_type = cols(Treated= col_logical()))%>%
    rename("depth_below_plot_cm"="Depth below plot base (cm)", 
         "depth_below_ground_cm"="Depth below ground surface (cm)") %>%
    column_to_rownames(.,var="SampleID")%>%
    sample_data()

#make featureID a row and the samples w/ # of reads and phyloseq otu object
table <- read_qza("table.qza")$data
otu<-otu_table(table, taxa_are_rows = TRUE)
  
#calculate relative abundance 
relabunKTYA <- as_tibble(table, rownames="featureID")%>%
  mutate_at(.,vars(-featureID), .funs = ~./sum(.)) 

#pull in taxonomy table, filter out lower confidence variables (<.9) then only pull taxon and featureID
taxonomy <- read_qza("taxonomy.qza")$data %>%
  as_tibble() %>%
  select (-Confidence) %>%
  mutate(Taxon=str_replace_all(string=Taxon,pattern = "[dpcofg]__", replacement = ""))%>%
  separate(Taxon,c("Domain","Phylum","Class","Order","Family","Genus"), sep = "; ")%>%
  column_to_rownames(.,var="Feature.ID")
#make the phyloseq version of the new taxonomy table. 
tax<- as.matrix(taxonomy)%>%
  tax_table()
taxonomy<- rownames_to_column(taxonomy, var= "Feature.ID")
#Download rootedtree
rootedtree <- read_qza("rooted-tree.qza")$data%>%
  phy_tree()


#make phyloseq obj
phylo<- phyloseq(otu,tax,metadata, rootedtree)



#join taxonomy and table so that taxonomy comes first. mutate and split taxon down to genus lecvel
#linked_KTYA <- inner_join(taxonomy_KTYA,table_KTYA, by = c("Feature.ID"="featureID") )%>%
 # mutate(Taxon=str_replace_all(string=Taxon,pattern = "[dpcofg]__", replacement = ""))%>%
  #select(-Feature.ID)%>%
  #separate(Taxon,c("Domain","Phylum","Class","Order","Family","Genus"), sep = ";")


#join relative abundance same as table
#linked_RA_KTYA <- inner_join(taxonomy_KTYA,relabunKTYA, by = c("Feature.ID"="featureID") )%>%
  #mutate(Taxon=str_replace_all(string=Taxon,pattern = "[dpcofrmg]__", replacement = ""))%>%
  #select(-Feature.ID)




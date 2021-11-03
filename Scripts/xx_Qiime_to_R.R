#This script pulls in files made in qiime and sets them to be used by R
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(here)

#renaming columns with problematic names
metadata<-read_tsv("Environment/metadata.txt", col_type = cols(Treated= col_logical(),
                              Timing = col_factor(levels = c("Before","After")), 
                              Treatment = col_factor(levels = c("NS","MX","WC","BC")),
                              Location = col_factor(levels = c("KTYA","KTR","HSP"))))
saveRDS(metadata, "Objects/metadata.rds")
   
metadata_phyloseq<- column_to_rownames(metadata,var="SampleID")%>% 
  sample_data()

#make featureID a row and the samples w/ # of reads and phyloseq otu object
table <- read_qza("Qiime/table.qza")$data 
otu<-otu_table(table, taxa_are_rows = TRUE)
 a<-  sample_names(otu)%>%
   as_tibble()
#calculate relative abundance 
relabun <- as_tibble(table, rownames="featureID")%>%
  mutate_at(.,vars(-featureID), .funs = ~./sum(.)) 

#pull in taxonomy table, filter out lower confidence variables (<.9) then only pull taxon and featureID
taxonomy <- read_qza("Qiime/taxonomy.qza")$data %>%
  as_tibble() %>%
  select (-Confidence) %>%
  mutate(Taxon=str_replace_all(string=Taxon,pattern = "[dpcofgs]__", replacement = ""))%>%
  separate(Taxon,c("Domain","Phylum","Class","Order","Family","Genus", "Species"), sep = "; ")%>%
  column_to_rownames(.,var="Feature.ID")
saveRDS(taxonomy, "Objects/taxonomy.rds")
#make the phyloseq version of the new taxonomy table. 
tax<- as.matrix(taxonomy)%>%
  tax_table()
taxonomy<- rownames_to_column(taxonomy, var= "feature.ID")
saveRDS(taxonomy, "Objects/taxonomy.rds")
#Download rootedtree
rootedtree <- read_qza("Qiime/rooted-tree.qza")$data%>%
  phy_tree()


#make phyloseq obj
phylo<- phyloseq(otu,tax,metadata_phyloseq, rootedtree)
saveRDS(phylo, file = "Objects/phlyoseq.rds")

table <- as_tibble(table, rownames= "feature.ID")


table <- inner_join(taxonomy, table, by ="feature.ID")%>%
   select(-feature.ID,-Domain)


#join taxonomy and table so that taxonomy comes first. mutate and split taxon down to genus lecvel
#linked_KTYA <- inner_join(taxonomy_KTYA,table_KTYA, by = c("Feature.ID"="featureID") )%>%
 # mutate(Taxon=str_replace_all(string=Taxon,pattern = "[dpcofg]__", replacement = ""))%>%
  #select(-Feature.ID)%>%
  #separate(Taxon,c("Domain","Phylum","Class","Order","Family","Genus"), sep = ";")


#join relative abundance same as table
#linked_RA_KTYA <- inner_join(taxonomy_KTYA,relabunKTYA, by = c("Feature.ID"="featureID") )%>%
  #mutate(Taxon=str_replace_all(string=Taxon,pattern = "[dpcofrmg]__", replacement = ""))%>%
  #select(-Feature.ID)




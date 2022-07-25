library(pacman)

p_load(phyloseq,here,tidyverse, ggplot2, forcats, ggpubr, vegan, ggvegan, ggtext, ANCOMBC)


#reading in Phylum data
metadata <- readRDS(here("Objects/metadata.rds"))%>%
  subset(., Treatment %in% c("NS", "WC"))
metadata$Location<- factor(metadata$Location, levels = c("HSP", "KTR", "KTYA"), ordered = T)

#metadata$Timing<- factor(metadata$Timing, levels = c("Before", "After"))
#metadata$Treatment<- factor(metadata$Treatment, levels = c("NS", "MX","WC", "BC"))
phylo<- readRDS(here("Objects/genera.rds"))
phylo<-tax_glom(phylo,taxrank="Phylum")
#removing BC from samples
phylo<-subset_samples(phylo,Treatment %in% c("NS", "WC"))
#**saved RDS as phylum**

#subset to only look at before sample
#phylo<-subset_samples(phylo, Timing == "Before")
#add a group column that combinies timing wiht treatment
sample_data(phylo)$group <-factor(paste0(sample_data(phylo)$Timing,":",sample_data(phylo)$Treatment))
palette<-readRDS(here("Objects/palette.RDS"))
n_palette<-colorRampPalette(palette)(11)
#phylo_ra<-transform_sample_counts(phylo, function(x) x / sum(x))
#putting everything but top 5
top10<- names(sort(taxa_sums(phylo),TRUE)[1:10])

#phylo_top10<- prune_taxa(top10, phylo)

dat<-psmelt(phylo)%>%
  mutate(., Phylum= ifelse(OTU %in% top10, Phylum, "Other"))%>%
  group_by(group, Location) %>%
  mutate(gAbun= sum(Abundance)) %>%
  group_by(Phylum, .add=TRUE) %>%
  mutate(RA=Abundance/gAbun)
dat$Location<- factor(dat$Location, levels = c("HSP", "KTR", "KTYA"), ordered = T)
names<- dat%>%
  group_by(group, Timing, Location, Phylum)%>%
  mutate(total=sum(RA))%>%
  select(Location, Timing, Phylum, total)%>%
  distinct()

table<-as(otu_table(phylo), "matrix")

shannon<-diversity(table, index = "shannon")%>%
  as_tibble(rownames =  "sample")%>%
  inner_join(., metadata, by=c("sample"="SampleID"))

palette1<- c("#89C0DB" , "#7D9FAE" , "#655F54" , "#593F28" , "#476129" , "#666E4B" ,  "#B8B89F" ,"#6C855C" ,
            "#8AA67E" , "#A8937E" , "grey")


diversity<-ggplot(shannon)+
  geom_boxplot(aes(x=Timing, y=value, col=Timing))+
  facet_grid(cols = vars(Location))+
  scale_color_manual(values=c(palette[2], palette[1]))+
  theme_pubr(legend = "none")+
  theme(axis.text.x = element_text(angle=45, vjust=.4))+
  labs(y="Shannon Weaver Index", x="")+
  stat_compare_means(aes(x=Timing, y=value, label = ..p.signif..), method= "t.test", label.x = 1.5, label.y = 2.2,
              hide.ns = T)
#ggsave("Figures/shannon.png",width = 6, height = 6)
dat$Phylum<-with(dat, ifelse(Phylum == "Other", "Other", 
                            paste0("*", Phylum, "*")))
RA<-ggplot(dat, aes(x=Treatment, y=RA, fill=Phylum))+
  facet_grid(cols= vars(Location),rows= vars(Timing), scales = "free_x", space = "free_x")+
  geom_col(aes( fill=Phylum), position="stack")+
  scale_fill_manual(values=palette1)+
  theme_pubr(legend = "right")+
  #geom_text(aes(label=total))+
  labs(y="Relative Abundance",x="Treatment" )+
  theme(panel.grid.major.x = element_blank(),legend.text = element_markdown())

#pco1 <- wcmdscale(vegdist(t(table)), eig = TRUE,add = "lingoes")
#  inner_join(.,metadata, by="SampleID")
#ggplot(pcoa)+
 # geom_point(aes(x=Dim2, y=Dim3, col=Location, shape=Timing))+
  #scale_color_manual(values=palette)+
  #theme_pubr(legend = "right")
set.seed(3618)
ndms<-metaMDS(table, distance = "bray", try=40, autotransform = F)

f_ndms<- fortify(ndms)%>%
  subset(., Score=="sites")%>%
 left_join(., metadata, by=c("Label"="SampleID"))
anova<- adonis2(table~Location, data = metadata, 
                permutations = 10000 )
anova1<-adonis2(table~Location+Timing, data=metadata, permutations = 1000)
anova1
#doesn't significantly separate by timing

anovaT<- adonis2(table~Timing, data = metadata)
anovaT
f_ndms$Location<-factor(f_ndms$Location, levels = c("HSP", "KTR", "KTYA"), ordered = T)
bray<-ggplot(f_ndms)+
  geom_point(aes(x=NMDS1,y=NMDS2, col=Location, shape=Timing))+
  theme_pubr(legend = "right")+
  scale_color_manual(values=palette)+
    scale_fill_manual(values=palette)+
   scale_shape_manual(values = c(20,4, 25))+ 
  geom_abline(intercept = 0,slope = 0,linetype="dashed", size=0.8)+ 
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.8)+
  stat_ellipse( aes(x=NMDS1,y=NMDS2, col=Location, fill=Location, linetype=Location),
               geom="polygon", type="t", alpha= 0.2, show.legend = T, level = .95)+
    labs(caption = "PERMANOVA(Location) <0.0001\nPERMANOVA(Timing)=0.052\nStress=0.07")

bray

ggarrange(RA, ggarrange(diversity, bray, labels=c("B","C"), ncol=2),  
        labels = c("A", "B"),
         ncol = 1, nrow = 2,heights =c(1.5,1))
#have to manually save bc of lines
# testing if proteobacteria increas is significant
taxa<-as(tax_table(phylo), "matrix")%>%
  as_tibble(rownames = "Label")
hsp_phylo<-subset_samples(phylo, Location=="HSP")
KTR_phylo<-subset_samples(phylo, Location=="KTR")
KTYA_phylo<-subset_samples(phylo_top10, Location=="KTYA")
hspout = ancombc(phyloseq= hsp_phylo, formula = "Timing"  , p_adj_method = "holm", zero_cut=.9,
              lib_cut = 0, struc_zero=FALSE, neg_lb = FALSE,  alpha=0.05) 
KTRout = ancombc(phyloseq= KTR_phylo, formula = "Timing"  , p_adj_method = "holm", zero_cut=.9,
                 lib_cut = 0, struc_zero=FALSE, neg_lb = FALSE,  global=F,alpha=0.05) 
KTYAout = ancombc(phyloseq= KTYA_phylo, formula = "Timing" , p_adj_method = "holm", zero_cut=.9,
                 lib_cut = 50, struc_zero=FALSE, neg_lb = FALSE,  global=F, alpha=0.05) 

hspout <- data.frame(
  Species = row.names(hspout$res$beta),
  coef = unlist(hspout$res$beta),
  se = unlist(hspout$res$se),
  W = unlist(hspout$res$W),
  p_val = unlist(hspout$res$p_val),
  q_val = unlist(hspout$res$q_val),
  diff_abn = unlist(hspout$res$diff_abn))
#filter only  significant genera and change from natural log to 2-log
hspout<- subset(hspout, diff_abn==TRUE)%>%
  mutate(., log2= coef * log2(exp(1)))%>%
  mutate(., selog2= se *log2(exp(1)))%>%
  left_join(., taxa, by= c("Species" = "Label"))

KTRout <- data.frame(
  Species = row.names(KTRout$res$beta),
  coef = unlist(KTRout$res$beta),
  se = unlist(KTRout$res$se),
  W = unlist(KTRout$res$W),
  p_val = unlist(KTRout$res$p_val),
  q_val = unlist(KTRout$res$q_val),
  diff_abn = unlist(KTRout$res$diff_abn))
#filter only  significant genera and change from natural log to 2-log
KTRout<- KTRout%>%
  mutate(., log2= coef * log2(exp(1)))%>%
  mutate(., selog2= se *log2(exp(1)))%>%
  left_join(., taxa, by= c("Species" = "Label"))
KTYAout <- data.frame(
  Species = row.names(KTYAout$res$beta),
  coef = unlist(KTYAout$res$beta),
  se = unlist(KTYAout$res$se),
  W = unlist(KTYAout$res$W),
  p_val = unlist(KTYAout$res$p_val),
  q_val = unlist(KTYAout$res$q_val),
  diff_abn = unlist(KTYAout$res$diff_abn))
#filter only  significant genera and change from natural log to 2-log
KTYAout<- subset(KTYAout, diff_abn==TRUE)%>%
  mutate(., log2= coef * log2(exp(1)))%>%
  mutate(., selog2= se *log2(exp(1)))%>%
  left_join(., taxa, by= c("Species" = "Label"))


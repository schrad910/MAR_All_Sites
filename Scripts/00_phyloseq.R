library(phyloseq)
library(Biostrings)
library(phangorn)
library(DECIPHER)
#load genus- glommed phyloseq object
ps <- readRDS("~/Documents/SaltikovLab/Sequencing/MAR_All_Sites/Objects/genera.rds")

#making a tree and adding it to the phyloseq object
aligngenus<- AlignSeqs(refseq(genus), anchor=NA)
gphangAlign <- phyDat(as(aligngenus, "matrix"), type="DNA")
dm <- dist.ml(gphangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=gphangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
phy_tree(ps)<-phy_tree(fitGTR$tree)
#for unifrac adn faiths pd need a root
set.seed(318)
phy_tree(ps)<-root(phy_tree(ps),sample(taxa_names(ps),1),resolve.root = TRUE)
is.rooted(phy_tree(ps))

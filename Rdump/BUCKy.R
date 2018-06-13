# SCRIPT IS DESIGNED TO VISUALISE THE OUTPUT OF BUCKy

library(ape)
library(ggtree)


#Load all the trees from Bucky analysis

PT <- read.tree("/Users/sa01fd/Downloads/RE__bucky_out/population_tree.nwk")
PTBL <- read.tree("/Users/sa01fd/Downloads/RE__bucky_out/population_tree_branch_lengths.nwk")
PCT <- read.tree("/Users/sa01fd/Downloads/RE__bucky_out/primary_concordance_tree.nwk")
PCTCF <- read.tree("/Users/sa01fd/Downloads/RE__bucky_out/primary_concordance_tree_concordance_factors.nwk")
tip_labels<-read.csv("/Users/sa01fd/Downloads/RE__bucky_out/tip_labels.csv",header=F,row.names=1,check.names=FALSE)

#Change the labels of the trees
PT$tip.label<-as.character(tip_labels[PT$tip.label,])
PTBL$tip.label<-as.character(tip_labels[PTBL$tip.label,])
PCT$tip.label<-as.character(tip_labels[PCT$tip.label,])
PCTCF$tip.label<-as.character(tip_labels[PCTCF$tip.label,])



#Visualise the trees using ggtree

p1 = ggtree(PT) + geom_text2(aes(subset=!isTip, label=node), hjust=1.2,vjust=-0.6,size=2) + geom_tiplab(size=2)+geom_treescale(fontsize=2)+xlim(0,50)+ggtitle("Population Tree")
p2 = ggtree(PTBL) + geom_text2(aes(subset=!isTip, label=node), hjust=1.2,vjust=-0.6,size=2) + geom_tiplab(size=2)+geom_treescale(fontsize=2)+xlim(0,50)+ggtitle("Population Tree, With Branch Lengths In Estimated Coalescent Units")
p3 = ggtree(PCT) + geom_text2(aes(subset=!isTip, label=node), hjust=1.2,vjust=-0.6,size=2) + geom_tiplab(size=2)+geom_treescale(fontsize=2)+xlim(0,50)+ggtitle("Primary Concordance Tree Topology")
p4 = ggtree(PCTCF) + geom_text2(aes(subset=!isTip, label=node), hjust=1.2,vjust=-0.6,size=2) + geom_tiplab(size=2)+geom_treescale(fontsize=2)+xlim(0,10)+ggtitle("Primary Concordance Tree with Sample Concordance Factors")

#individual graphics
print(p1)
print(p2)
print(p3)
print(p4)

#combined graphic
multiplot(p1,p2,p3,p4,cols=2)



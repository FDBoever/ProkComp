# 	Frederik De Boever, 01-Jun-2018
#
# 	MultiTreeStat.R
#	MiltiTreeStat.R wraps basic visualisation and statistics of phylogenies (or other cluster obtjects) in R, using a variaty of establisched packages
#	
#	This Script is designed to 
#	1) 	visual inspection of trees at different K
#	2) 	visual inspection of tree comparisons using tanglegrams
#	3) 	calculate Fowlkes-Mallows index between two trees; at different K
#	4)	Correlation among a list of trees using, cophenetic, baker, common nodes and F-index
#	5)	Calculate Baker's Gamma correlation coefficient for two trees (also known as Goodman-Kruskal-gamma index).
#
#	Reference:
#	Fowlkes, E. B.; Mallows, C. L. (1 September 1983). "A Method for Comparing Two Hierarchical Clusterings". Journal of the American Statistical Association 78 (383): 553.
#	Baker, F. B., Stability of Two Hierarchical Grouping Techniques Case 1: Sensitivity to Data Errors. Journal of the American Statistical Association, 69(346), 440 (1974).
#########################################
#			INSTALL PACKAGES
#########################################

#source("https://bioconductor.org/biocLite.R")
#biocLite("treeio")
#devtools::install_github("GuangchuangYu/treeio")

#########################################
#			LOAD PACKAGES
#########################################
library(treeio)
library(ggtree)
library(dendextend)
library(corrplot)
library(phytools)
library(reshape)
library(corrplot)

#########################################
#	optional RAxML specific import
#########################################
#raxml_file <- system.file("extdata/RAxML", "~/Downloads/ds-2/RAxML_bipartitions.Campbell_AA", package="treeio")
#raxml_Campbell_AA <- read.raxml(raxml_file)
#raxml_file <- system.file("extdata/RAxML", "~/Downloads/ds-2/RAxML_bipartitions.Ribosomal_AA", package="treeio")
#raxml_Ribosomal_Proteins_AA <- read.raxml(raxml_file)
#raxml_file <- system.file("extdata/RAxML", "~/Downloads/ds-2/RAxML_bipartitions.Rinke_AA", package="treeio")
#raxml_Rinke_AA <- read.raxml(raxml_file)
#raxml_file <- system.file("extdata/RAxML", "~/Downloads/ds-2/SCO_682_RAxML_bipartitions_autosubst_b100", package="treeio")
#raxml_ANVIO_SCO_AA <- read.raxml(raxml_file)

# optional meta_data input
phyloMETA2 = data.frame(cbind("set"=c("SCO","Campbell","Rinke","RP"),"columns"=c(123465,50908,32758,6844),"informative"=c(48794,13114,7909,1393),"constant"=c(62732,33799,21910,5164)))
phyloMETA$value=as.numeric(as.character(phyloMETA$value))

phyloMETA2$columns =as.numeric(as.character(phyloMETA2$columns))
phyloMETA2$informative =as.numeric(as.character(phyloMETA2$informative))
phyloMETA2$constant =as.numeric(as.character(phyloMETA2$constant))

phyloMETA = melt(phyloMETA2,id="set")


ggplot(phyloMETA,aes(reorder(set,value),value,fill= variable))+geom_bar(stat='identity',position=position_dodge())+scale_fill_manual(values=c('black','darkgray','lightgray'))+theme_classic()+facet_wrap(~variable)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_y_continuous(expand = c(0,0))+theme(strip.placement = "outside",strip.background = element_blank(),strip.text = element_text(face = "bold"))+xlab('marker sets')

ggplot(phyloMETA,aes(reorder(set,value),value,fill= variable))+geom_bar(stat='identity',position=position_dodge())+scale_fill_manual(values=c('black','darkgray','lightgray'))+theme_classic()+facet_wrap(~variable,scale='free_y')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_y_continuous(expand = c(0,0))+theme(strip.placement = "outside",strip.background = element_blank(),strip.text = element_text(face = "bold"))+xlab('marker sets')

#########################################
#	COLORS
#########################################

colors = c(
rgb(57,129,29, maxColorValue=255),#adhaerens
rgb(134,218,129, maxColorValue=255),#alg
rgb(166,216,212, maxColorValue=255),#ant
rgb(242,169,104, maxColorValue=255),#ex

rgb(242,236,112, maxColorValue=255), #hydrocarb
rgb(227,143,221, maxColorValue=255),#lipo

rgb(137,137,137, maxColorValue=255),#other
rgb(118,174,207, maxColorValue=255),#psych
rgb(131,115,27, maxColorValue=255),#sedim
rgb(179,77,34, maxColorValue=255)#vinifirmus
)

Habitat_colors = c('gray',"#E41A1C","#F37912","#FFD422","#43997A","#658E67","#5D6795","#A35390","#B85F49")


#########################################
#	tree import and dendogram transformation
#########################################

#raxml_ANVIO_SCO_AA = read.tree("~/Downloads/ds-2/SCO_682_RAxML_bipartitions_autosubst_b100")
raxml_ANVIO_SCO_AA = read.tree("~/DATA/MarinobacterGenomics/2018_ProkComp/trees/SCO_kde.fas.treefile")
midRoot_raxml_ANVIO_SCO_AA = ladderize(midpoint.root(raxml_ANVIO_SCO_AA))
dendSCO = chronos(midRoot_raxml_ANVIO_SCO_AA)
dendSCO = as.dendrogram(dendSCO)

raxml_Ribosomal_Proteins_AA = read.tree("~/DATA/MarinobacterGenomics/2018_ProkComp/trees/concatenated_Ribosomal_proteins_AA.fa.bmge.treefile")
#raxml_Ribosomal_Proteins_AA = read.tree("~/Downloads/ds-2/RAxML_bipartitions.Ribosomal_AA")
midRoot_raxml_Ribosomal_Proteins_AA = ladderize(midpoint.root(raxml_Ribosomal_Proteins_AA))
dendRP = chronos(midRoot_raxml_Ribosomal_Proteins_AA)
dendRP = as.dendrogram(dendRP)

#raxml_Campbell_AA = read.tree("~/Downloads/ds-2/RAxML_bipartitions.Campbell_AA")
raxml_Campbell_AA = read.tree("~/DATA/MarinobacterGenomics/2018_ProkComp/trees/concatenated_Campbell_AA.fa.bmge.treefile")
midRoot_raxml_Campbell_AA = ladderize(midpoint.root(raxml_Campbell_AA))
dendCampbell = chronos(midRoot_raxml_Campbell_AA)
dendCampbell = as.dendrogram(dendCampbell)

#raxml_Rinke_AA = read.tree("~/Downloads/ds-2/RAxML_bipartitions.Rinke_AA")
raxml_Rinke_AA = read.tree("~/DATA/MarinobacterGenomics/2018_ProkComp/trees/concatenated_Rinke_AA.fa.bmge.treefile")
midRoot_raxml_Rinke_AA = ladderize(midpoint.root(raxml_Rinke_AA))
dendRinke = chronos(midRoot_raxml_Rinke_AA)
dendRinke = as.dendrogram(dendRinke)

#########################################
#	tree import and dendogram transformation
#########################################
ANVIO_cat = read.table('~/DATA/MarinobacterGenomics/2018_ProkComp/ANVIO_CAT.txt',header=TRUE,sep="\t")
ANVIO_cat = ANVIO_cat[grepl(paste(midRoot_raxml_ANVIO_SCO_AA$tip.label,collapse='|'), ANVIO_cat$name),]
rownames(ANVIO_cat) = ANVIO_cat$name
ANVIO_cat  = ANVIO_cat[midRoot_raxml_ANVIO_SCO_AA$tip.label,]
ANVIO_cat[ANVIO_cat$SS=='sediment','SS'] = 'Sediment'



t1 = ggtree(midRoot_raxml_ANVIO_SCO_AA) + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_point2(aes(subset=!isTip), color='black', size=2)+ geom_tiplab(size=2) + xlim(0,2)

t1 = ggtree(midRoot_raxml_ANVIO_SCO_AA) + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(size=2) + xlim(0,1)
t2 = ggtree(midRoot_raxml_Campbell_AA) + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(size=2)+ xlim(0,1)
t3 = ggtree(midRoot_raxml_Rinke_AA) + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) +  geom_tiplab(size=2)+ xlim(0,1)
t4 = ggtree(midRoot_raxml_Ribosomal_Proteins_AA) + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(size=2)+ xlim(0,0.5)

multiplot(t1,t2,t3,t4,cols=4)

t1 =ggtree(midRoot_raxml_ANVIO_SCO_AA) %<+% ANVIO_cat +  geom_tiplab(aes(color = group2)) + geom_point(aes(color = group2)) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2)+ xlim(0,1)
t2 = ggtree(midRoot_raxml_Campbell_AA) %<+% ANVIO_cat +  geom_tiplab(aes(color = group2)) + geom_point(aes(color = group2)) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2)+ xlim(0,1)
t3 = ggtree(midRoot_raxml_Rinke_AA) %<+% ANVIO_cat +  geom_tiplab(aes(color = group2)) + geom_point(aes(color = group2)) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2)+ xlim(0,1)
t4 =ggtree(midRoot_raxml_Ribosomal_Proteins_AA) %<+% ANVIO_cat +  geom_tiplab(aes(color = group2)) + geom_point(aes(color = group2)) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2)+ xlim(0,1)

multiplot(t1,t2,t3,t4,cols=4)


t1 = ggtree(midRoot_raxml_ANVIO_SCO_AA,branch.length = 'none' ) %<+% ANVIO_cat + geom_point(aes(color = group2),size=3) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2)
t2 = ggtree(midRoot_raxml_Campbell_AA,branch.length = 'none' ) %<+% ANVIO_cat + geom_point(aes(color = group2),size=3) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2)
t3 = ggtree(midRoot_raxml_Rinke_AA,branch.length = 'none' ) %<+% ANVIO_cat + geom_point(aes(color = group2),size=3) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2)
t4 = ggtree(midRoot_raxml_Ribosomal_Proteins_AA,branch.length = 'none' ) %<+% ANVIO_cat + geom_point(aes(color = group2),size=3) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2)


#color coded based on habitat
t1.1 = ggtree(midRoot_raxml_ANVIO_SCO_AA,branch.length = 'none' ) %<+% ANVIO_cat + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) +geom_tippoint(aes(color= SS), alpha=1,size=3) +  scale_color_manual(values= Habitat_colors)#+theme(legend.position='right')+geom_tiplab()
t2.1 = ggtree(midRoot_raxml_Campbell_AA,branch.length = 'none' ) %<+% ANVIO_cat + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) +geom_tippoint(aes(color= SS), alpha=1,size=3) +  scale_color_manual(values= Habitat_colors)#+theme(legend.position='right')+geom_tiplab()
t3.1 = ggtree(midRoot_raxml_Rinke_AA,branch.length = 'none' ) %<+% ANVIO_cat + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) +geom_tippoint(aes(color= SS), alpha=1,size=3) +  scale_color_manual(values= Habitat_colors)#+theme(legend.position='right')+geom_tiplab()
t4.1 = ggtree(midRoot_raxml_Ribosomal_Proteins_AA,branch.length = 'none' ) %<+% ANVIO_cat + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) +geom_tippoint(aes(color= SS), alpha=1,size=3) +  scale_color_manual(values=Habitat_colors)#+theme(legend.position='right')+geom_tiplab()
multiplot(t1.1,t2.1,t3.1,t4.1,cols=4)


t1 = ggtree(midRoot_raxml_ANVIO_SCO_AA,branch.length = 'none',layout="circular" ) %<+% ANVIO_cat + geom_point(aes(color = group2),size=3) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(offset=1.2, size=2,aes(angle=angle))+geom_tippoint(aes(color= group2), alpha=1,size=2)

t12 = ggtree(midRoot_raxml_ANVIO_SCO_AA,branch.length = 'none' ) %<+% ANVIO_cat +  scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(offset=1.2, size=2,aes(angle=angle))+geom_tippoint(aes(color= group2), alpha=1,size=3)
t12 = open_tree(t12, angle=180)
t13 = ggtree(midRoot_raxml_Campbell_AA,branch.length = 'none' ) %<+% ANVIO_cat +  scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(offset=1.2, size=2,aes(angle=angle))+geom_tippoint(aes(color= group2), alpha=1,size=3)
t13 = open_tree(t13, angle=180)
t14 = ggtree(midRoot_raxml_Rinke_AA,branch.length = 'none' ) %<+% ANVIO_cat +  scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(offset=1.2, size=2,aes(angle=angle))+geom_tippoint(aes(color= group2), alpha=1,size=3)
t14 = open_tree(t14, angle=180)
t15 = ggtree(midRoot_raxml_Ribosomal_Proteins_AA,branch.length = 'none' ) %<+% ANVIO_cat +  scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(offset=1.2, size=2,aes(angle=angle))+geom_tippoint(aes(color= group2), alpha=1,size=3)
t15 = open_tree(t15, angle=180)

multiplot(t12,t13,t14,t15)

t12 = ggtree(midRoot_raxml_ANVIO_SCO_AA,branch.length = 'none' ) %<+% ANVIO_cat +  scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(offset=1.2, size=2,aes(angle=angle))+geom_tippoint(aes(color= group2), alpha=1,size=3)
t12 = open_tree(t12, angle=180)+theme(legend.position='bottom')
t12.1 = ggtree(midRoot_raxml_ANVIO_SCO_AA,branch.length = 'none' ) %<+% ANVIO_cat +  scale_colour_manual(values= Habitat_colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(offset=1.2, size=2,aes(angle=angle))+geom_tippoint(aes(color= SS), alpha=1,size=3)
t12.1 = open_tree(t12.1, angle=180)+theme(legend.position='bottom')



ANVIO_cat$GC = traitGC[rownames(ANVIO_cat)]
ANVIO_cat$size = traitSize[rownames(ANVIO_cat)]

t12.2 = ggtree(midRoot_raxml_ANVIO_SCO_AA,branch.length = 'none' ) %<+% ANVIO_cat + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(offset=1.2, size=2,aes(angle=angle))+geom_tippoint(aes(color= GC), alpha=1,size=3)+scale_colour_gradient(low = "white", high = "darkslategray")
t12.2 = open_tree(t12.2, angle=180)+theme(legend.position='bottom')

t12.3 = ggtree(midRoot_raxml_ANVIO_SCO_AA,branch.length = 'none' ) %<+% ANVIO_cat + geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(offset=1.2, size=2,aes(angle=angle))+geom_tippoint(aes(color= size), alpha=1,size=3)+scale_colour_gradient(low = "white", high = "burlywood4")
t12.3 = open_tree(t12.3, angle=180)+theme(legend.position='bottom')


t12.1 = ggtree(midRoot_raxml_ANVIO_SCO_AA,branch.length = 'none' ) %<+% ANVIO_cat +  scale_colour_manual(values= Habitat_colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2) + geom_tiplab(offset=1.2, size=2,aes(angle=angle))+geom_tippoint(aes(color= SS), alpha=1,size=3) 
t12.1 <-  gheatmap(t12.1, data.frame(traitGC), width=0.05, low="white", high="darkslategray", colnames = FALSE)
t12.1 = open_tree(t12.1, angle=180)+theme(legend.position='bottom')




#########################################
#	make a list of dendrograms
#########################################

dends=dendlist("SCO" = dendSCO, "Campbell" = dendCampbell, "Rinke" = dendRinke, "RP" = dendRP)
dends2=dendlist("SCO" = dendSCO, "Campbell" = dendCampbell)

#########################################
#	Bk plots
#########################################
#Reference:
#Fowlkes, E. B.; Mallows, C. L. (1 September 1983). "A Method for Comparing Two Hierarchical Clusterings". Journal of the American Statistical Association 78 (383): 553.#Fowlkes-Mallows index (see references) is an external evaluation method that is used to determine the similarity between two clusterings (clusters obtained after a clustering algorithm). This measure of similarity could be either between two hierarchical clusterings or a clustering and a benchmark classification. A higher the value for the Fowlkes-Mallows index indicates a greater similarity between the clusters and the benchmark classifications.
#The default Bk plot comes with a line with dots (type "b") of the Bk values. Also with a fragmented (lty=2) line (of the same color) of the expected Bk line under H0, And a solid red line of the upper critical Bk values for rejection
#Bk Plot - Ploting The Fowlkes-Mallows Index Of Two Dendrogram For Various K's
#Bk is the calculation of Fowlkes-Mallows index for a series of k cuts for two dendrograms. A Bk plot is simply a scatter plot of Bk versus k. This plot helps in identifiying the similarity between two dendrograms in different levels of k (number of clusters).

#without permutation
par(mfrow = c(3,3))
Bk_plot(dendSCO, dendCampbell)
Bk_plot(dendSCO, dendRinke)
Bk_plot(dendSCO, dendRP)
Bk_plot(dendCampbell, dendSCO)
Bk_plot(dendCampbell, dendRinke)
Bk_plot(dendCampbell, dendRP)
Bk_plot(dendRinke, dendSCO)
Bk_plot(dendRinke, dendCampbell)
Bk_plot(dendRinke, dendRP)

#with permutation
par(mfrow = c(3,3))
Bk_plot(dendSCO, dendCampbell,rejection_line_permutation=TRUE,  R= 100)
Bk_plot(dendSCO, dendRinke,rejection_line_permutation=TRUE,  R= 100)
Bk_plot(dendSCO, dendRP,rejection_line_permutation=TRUE,  R= 100)
Bk_plot(dendCampbell, dendSCO,rejection_line_permutation=TRUE,  R= 100)
Bk_plot(dendCampbell, dendRinke,rejection_line_permutation=TRUE,  R= 100)
Bk_plot(dendCampbell, dendRP,rejection_line_permutation=TRUE,  R= 100)
Bk_plot(dendRinke, dendSCO,rejection_line_permutation=TRUE,  R= 100)
Bk_plot(dendRinke, dendCampbell,rejection_line_permutation=TRUE,  R= 100)
Bk_plot(dendRinke, dendRP,rejection_line_permutation=TRUE,  R= 100)

#########################################
#	visual inspection of trees at different K
#########################################

par(mfrow = c(4,4),oma = c(5,4,0,0) + 0.1,mar = c(0,0,1,1) + 0.1))
for(o in c(2,5,12,20)) {
   for(i in 1:length(dends)) {
   dends[[i]] %>% set("branches_k_color", k=o) %>% set("labels_cex", 0.2)  %>% plot(axes = FALSE, horiz = TRUE)
   title(paste(names(dends)[i], " k=",o))
   }
}

#########################################
#	Correlation among a list of trees
#########################################

par(mfrow = c(1,1))
all.equal(dends)
all.equal(dends, use.edge.length = FALSE)
cor.dendlist(dends)
corrplot(cor.dendlist(dends), "pie", "lower")

par(mfrow = c(2,2))
corrplot(cor.dendlist(dends,method="cophenetic"), "pie", "lower",main='cophenetic')
corrplot(cor.dendlist(dends,method="baker"), "pie", "lower",main='baker')
corrplot(cor.dendlist(dends,,method="common_nodes"), "pie", "lower",main='common_nodes')
corrplot(cor.dendlist(dends,method="FM_index",k=10), "pie", "lower",main="FM_index, k=10")


#########################################
#	comparison by tanglegrams
#########################################

par(mfrow = c(1,3),oma = c(5,4,0,0) + 0.1,mar = c(0,0,1,1) + 0.1))
for(i in 1:length(dends)-1) {
   dends2 =  dendlist(dends[[1]], dends[[i+1]])   
   dends2 %>%plot(common_subtrees_color_branches = TRUE, main = paste(names(dends)[1], names(dends)[i+1]))
   
}

#########################################
#	common leaves
#########################################

common_leaves = cbind()
for(i in 1:length(dends)-1) {
   	dends2 =  dendlist(dends[[1]], dends[[i+1]])   
   	d_common <- dends2 %>% prune_common_subtrees.dendlist 
	d_common  %>%  t%>%plot(common_subtrees_color_branches = TRUE, main = paste(names(dends)[1], names(dends)[i+1]))
	common_leaves = cbind(common_leaves,d_common %>% nleaves)  
}

#########################################
#	Baker's Gamma Correlation Coefficient
#########################################
#Calculate Baker's Gamma correlation coefficient for two trees (also known as Goodman-Kruskal-gamma index).
set.seed(23235)
the_cor2 <- cor_bakers_gamma(dendSCO, dendCampbell)
the_cor                    
the_cor2                
the_cor2
par(mfrow = c(1,1))

R <- 100
cor_bakers_gamma_results <- numeric(R)
dend_mixed <- dendSCO
for(i in 1:R) {
   dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)
   cor_bakers_gamma_results[i] <- cor_bakers_gamma(dendSCO, dend_mixed)
}
plot(density(cor_bakers_gamma_results),
     main = "Baker's gamma distribution under H0",
     xlim = c(-1,1))
abline(v = 0, lty = 2)
abline(v = the_cor, lty = 2, col = 2)
abline(v = the_cor2, lty = 2, col = 4)
legend("topleft", legend = c("cor", "cor2"), fill = c(2,4))
round(sum(the_cor2 < cor_bakers_gamma_results)/ R, 4)                    
title(sub = paste("One sided p-value:",
                  "cor =",  round(sum(the_cor < cor_bakers_gamma_results)/ R, 4),
                  " ; cor2 =",  round(sum(the_cor2 < cor_bakers_gamma_results)/ R, 4)
                  ))                
              
dend1 <- dendSCO
dend2 <- dendCampbell

set.seed(23801)

R <- 100
dend1_labels <- labels(dend1)
dend2_labels <- labels(dend2)
cor_bakers_gamma_results <- numeric(R)
for(i in 1:R) {
   sampled_labels <- sample(dend1_labels, replace = TRUE)
   # members needs to be fixed since it will be later used in nleaves
   dend_mixed1 <- sample.dendrogram(dend1, 
                                    dend_labels=dend1_labels,
                                    fix_members=TRUE,fix_order=TRUE,fix_midpoint=FALSE,
                                    replace = TRUE, sampled_labels=sampled_labels
                                      )
   dend_mixed2 <- sample.dendrogram(dend2, dend_labels=dend2_labels,
                                    fix_members=TRUE,fix_order=TRUE,fix_midpoint=FALSE,
                                    replace = TRUE, sampled_labels=sampled_labels
                                      )                                    
   cor_bakers_gamma_results[i] <- cor_bakers_gamma(dend_mixed1, dend_mixed2, warn = FALSE)
}
CI95 <- quantile(cor_bakers_gamma_results, probs=c(.025,.975))                  
par(mfrow = c(1,1))
plot(density(cor_bakers_gamma_results),
     main = "Baker's gamma bootstrap distribution",
     xlim = c(-1,1))
abline(v = CI95, lty = 2, col = 3)
abline(v = cor_bakers_gamma(dend1, dend2), lty = 2, col = 2)
legend("topleft", legend =c("95% CI", "Baker's Gamma Index"), fill = c(3,2))





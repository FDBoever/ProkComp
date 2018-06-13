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
phyloMETA2 = data.frame(cbind("set"=c("SCO","Campbell","Rinke","RP"),"columns"=c(177229,50858,32572,6841),"informative"=c(74378,14369,8783,1388),"constant"=c(87873,32952,21353,5065),"alignmentPattern"=c(78357,13565,8459,1698)))
phyloMETA = melt(phyloMETA2,id="set")
phyloMETA$value=as.numeric(as.character(phyloMETA$value))

phyloMETA2$columns =as.numeric(as.character(phyloMETA2$columns))
phyloMETA2$informative =as.numeric(as.character(phyloMETA2$informative))
phyloMETA2$constant =as.numeric(as.character(phyloMETA2$constant))
phyloMETA2$alignmentPattern =as.numeric(as.character(phyloMETA2$alignmentPattern))

#########################################
#	tree import and dendogram transformation
#########################################

raxml_ANVIO_SCO_AA = read.tree("~/Downloads/ds-2/SCO_682_RAxML_bipartitions_autosubst_b100")
midRoot_raxml_ANVIO_SCO_AA = ladderize(midpoint.root(raxml_ANVIO_SCO_AA))
dendSCO = chronos(midRoot_raxml_ANVIO_SCO_AA)
dendSCO = as.dendrogram(dendSCO)

raxml_Ribosomal_Proteins_AA = read.tree("~/Downloads/ds-2/RAxML_bipartitions.Ribosomal_AA")
midRoot_raxml_Ribosomal_Proteins_AA = ladderize(midpoint.root(raxml_Ribosomal_Proteins_AA))
dendRP = chronos(midRoot_raxml_Ribosomal_Proteins_AA)
dendRP = as.dendrogram(dendRP)

raxml_Campbell_AA = read.tree("~/Downloads/ds-2/RAxML_bipartitions.Campbell_AA")
midRoot_raxml_Campbell_AA = ladderize(midpoint.root(raxml_Campbell_AA))
dendCampbell = chronos(midRoot_raxml_Campbell_AA)
dendCampbell = as.dendrogram(dendCampbell)

raxml_Rinke_AA = read.tree("~/Downloads/ds-2/RAxML_bipartitions.Rinke_AA")
midRoot_raxml_Rinke_AA = ladderize(midpoint.root(raxml_Rinke_AA))
dendRinke = chronos(midRoot_raxml_Rinke_AA)
dendRinke = as.dendrogram(dendRinke)

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





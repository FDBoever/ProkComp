###############################################################################################
#
#	R script for PAN GENOME ANALYSIS
#
# 	(Script can use any presence/absence tables like Pfam domains)
#	in out study it was used as a downstream statistical pipeline for the ANVIO PAN GENOME ANALYSIS
#	it generates, pan-genome accululation curves, estimates closed/open pangenome, visualiation
#	weigthed distances, Jaccard and Manhattan, related Principal component analysis and network trees
#
###############################################################################################


#----------------------------------------------
#	installing packages
#----------------------------------------------

#install.packages("micropan",dependencies=TRUE)
#install.packages("phangorn",dependencies=TRUE)
#install.packages("ggplot2",dependencies=TRUE)
#install.packages("vegan",dependencies=TRUE)
#install.packages("devtools")
#library("devtools")
#install_github("ggbiplot", "vqv")
#install.packages("ade4")
#----------------------------------------------
#	Custum functions
#----------------------------------------------
generate_label_df <- function(HSD, flev){
 # Extract labels and factor levels from Tukey post-hoc 
 Tukey.levels <- HSD[[flev]][,4]
 Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
 plot.labels <- names(Tukey.labels[['Letters']])

 # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
 # upper quantile and label placement
    boxplot.df <- ddply(d, flev, function (x) max(fivenum(x$y)) + 0.2)

 # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
     stringsAsFactors = FALSE)

 # Merge it with the labels
   labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)

return(labels.df)
}


#----------------------------------------------
#	loading packages
#----------------------------------------------
library(micropan)
library(phangorn)
library(vegan)
library("ggbiplot")
library(ade4)
library(rgl)
library(ape)
library(phytools)
library(gplots)
library(ggplot2)
library(gridExtra)
library(rgl)
library(RColorBrewer)
library(reshape)
library(RColorBrewer)
library(dplyr)
library(seqinr)
library(viridis)
library(heatmap.plus)
devtools::install_github("chrislad/phenotypicForest")
library("multcompView")

#----------------------------------------------
#	Load PANGENOME DATA 
#----------------------------------------------

ANVIO <- read.table('~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/MARINOBACTER_protein_clusters_summary.txt',header=TRUE,sep="\t") #protein clusters
AnvioData <- read.table('~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PanGenome_matrix_pres_abs.txt',header=TRUE) #presence absence data
PfamData <- read.table('~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/Pfam_matrix_PA.txt',header=TRUE)
PfamData2 <- read.table('~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/Pfam_matrix.txt',header=TRUE)

ANVIO_cat <- read.table('~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/ANVIO_CAT.txt',header=TRUE,sep="\t")
CheckM = read.table("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/QC/qa2Marinobacter.txt",header=TRUE)
ANVIOcheck = read.table("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/QC/MARINOBACTER-samples-information.txt",header=TRUE)
BUSCOcheck = read.table("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/QC/Marinobacter_BUSCO_SUMMARY.txt",header = TRUE)
CompareCharacter = read.table("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/QC/CheckM_vs_others.txt",header = TRUE,sep="\t")



CheckM2 = subset(CheckM, Completeness > 98)
CheckM2 = CheckM2[order(CheckM2$Bin_Id),]

#Transpose matrices
AnvioData = t(AnvioData)
PfamData = t(PfamData)
PfamData2 = t(PfamData2)

#----------------------------------------------
#	Load Phylogenetic tree 
#----------------------------------------------
#load
RAxMLANVIO = read.tree("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/ANVIO_RAxML/RAxML_bipartitions.autosubst100")
#Root
RAxMLANVIORooted =  midpoint.root(RAxMLANVIO)
#increasing branch lengths
tree2 <- ladderize(RAxMLANVIORooted, right = FALSE)

#can’t have any branch lengths of zero or downstream commands will collapse those nodes…
tree2$edge.length[which(tree2$edge.length == 0)] <- 0.00001
rep_tree_um <- chronopl(tree2,lambda = 0.1,tol = 0)
rep_tree_d <- as.dendrogram(as.hclust.phylo(rep_tree_um))


cols <- RColorBrewer::brewer.pal(length(unique(ANVIO_cat$group)), name = "RdYlBu")
ANVIO_cat$color <- factor(ANVIO_cat$group, labels = cols)
X<-cbind(ANVIO_cat$lat, ANVIO_cat$lon)
rownames(X)=ANVIO_cat$name

obj<-phylo.to.map(tree2,X,ftype="i",fsize=0.7,
asp=1.2,split=c(0.5,0.5))


library(fossil)
GeoDist = earth.dist(X)


#----------------------------------------------
#	PC distribution
#----------------------------------------------

#count number of genomes in the dataset
number_of_genomes = length(colnames(AnvioData))

# number of genomes per PC and Pfam
sums = data.frame(colSums(AnvioData))
sumsPfam = data.frame(colSums(PfamData))
colnames(sums) = c('value')
colnames(sumsPfam) = c('value')

#Combine in one dataframe
groupedDistribution = cbind(rbind(sumsPfam,sums),type = c(rep("domain",nrow(sumsPfam)),rep("PC",nrow(sums))))
groupedDistribution = groupedDistribution[-c(which(groupedDistribution$value == 0)), ]

#Individual plots
ggplot(data=sums, aes(sums$value)) + 
  geom_histogram(bins = 60) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())
ggplot(data= sumsPfam, aes(sumsPfam $value)) + 
  geom_histogram(bins = 60) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())

o1 = ggplot(data= groupedDistribution, aes(x=value, fill=type)) + 
  geom_histogram(bins = 60,alpha=0.5, position="identity") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) +theme(axis.line = element_line(color="black", size = 0.5))
  
countTable <-table(groupedDistribution$value, groupedDistribution$type)
dfCount = data.frame(countTable)
countPC = subset(dfCount ,dfCount$Var2=="PC")
countDomain = subset(dfCount , dfCount$Var2 =="domain")
combinedCount = data.frame(cbind(PC=countPC$Freq,Domain=countDomain$Freq))

o2 = ggplot(data= combinedCount , aes(x=PC, y= Domain)) + geom_point() + scale_y_log10() +  scale_x_log10() + geom_smooth(method="lm",color="darkgreen") + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())+theme(axis.line = element_line(color="black", size = 0.5))

multiplot(o1,o2,cols=2)

#----------------------------------------------
#	extract accessory genes
#----------------------------------------------

AccesoryDF = AnvioData[,!colSums(AnvioData)==1 ]
AccesoryDF  = AccesoryDF[,!colSums(AccesoryDF)==60 ]
RareAccesoryDF  = AccesoryDF[,!colSums(AccesoryDF)<7]

#force row order so that it matches the order of leafs in rep_tree_d
clade_order <- order.dendrogram(rep_tree_d)
clade_name <- labels(rep_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(AccesoryDF))
combined_ordered_matrix <- AccesoryDF[new_order,]

heatmap.2(as.matrix(combined_ordered_matrix),trace="none",col =c('white',"black"),Rowv= rep_tree_d,labCol = FALSE,margins = c(2, 20),dendrogram="row")

combined_ordered_matrix <- RareAccesoryDF[new_order,]
heatmap.2(as.matrix(combined_ordered_matrix),trace="none",col =c('white',"black"),Rowv= rep_tree_d,labCol = FALSE,margins = c(2, 20),dendrogram="row")



#----------------------------------------------
#	distribution of unique genes
# sorted like phylogenetic tree! (plot against tree)
#----------------------------------------------

#ANVIO
TotalPCs = rowSums(AnvioData)
Core_count = rowSums(AnvioData[,colSums(AnvioData)==60 ])
only_singletons = AnvioData[,colSums(AnvioData)==1 ]
SummaryPanGenome = data.frame(rowSums(only_singletons))
colnames(SummaryPanGenome) = c('Unique')
SummaryPanGenome$name = rownames(SummaryPanGenome)
SummaryPanGenome$Core = Core_count 
SummaryPanGenome$Accessory = rowSums(AnvioData) - SummaryPanGenome$Core - SummaryPanGenome$Unique
SummaryPanGenome$total = SummaryPanGenome$Unique + SummaryPanGenome$Accessory + SummaryPanGenome$Core
SummaryPanGenome$conserved = SummaryPanGenome$Core / SummaryPanGenome$total 


#Pfam
TotalDomains = rowSums(PfamData)
Core_Pfam_count = rowSums(PfamData[,colSums(PfamData)==60 ])
Pfam_singletons = PfamData[,colSums(PfamData)==1 ]
SummaryPfamPanGenome = data.frame(rowSums(Pfam_singletons))
colnames(SummaryPfamPanGenome) = c('Unique')
SummaryPfamPanGenome$name = rownames(SummaryPfamPanGenome)
SummaryPfamPanGenome$Core = Core_Pfam_count 
SummaryPfamPanGenome$Accessory = rowSums(PfamData) - SummaryPfamPanGenome$Core - SummaryPfamPanGenome$Unique


k1 = ggplot(data= SummaryPanGenome, aes(y=SummaryPanGenome$Unique,x=SummaryPanGenome$name)) + geom_col() + xlab("Genome") + ylab("number of unique PCs") + coord_flip() + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())

k2 = ggplot(data= SummaryPfamPanGenome, aes(y= SummaryPfamPanGenome$Unique,x= SummaryPfamPanGenome$name)) + geom_col() + xlab("Genome") + ylab("number of unique Pfam domains") + coord_flip() + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())


multiplot(k1, k2,cols=2)


DF1 <- melt(SummaryPanGenome, id.var="name")
DF2 <- melt(SummaryPfamPanGenome, id.var="name")

p1 = ggplot(DF1, aes(x = name, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + xlab("Genome") + ylab("Protein Cluster") + coord_flip() + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())

p2 = ggplot(DF2, aes(x = name, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + xlab("Genome") + ylab("Pfam domains") + coord_flip() + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())

p3 = ggplot(subset(DF1,variable!='Core'), aes(x = variable, y = value)) + geom_boxplot() + xlab("pan-genome devision") + ylab("number of PCs") + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) + geom_hline(yintercept = unique(Core_count),linetype = 2) + theme_classic()

multiplot(p1,p2,cols=2)

Combined_distribution = cbind(DF1,DF2)
Combined_distribution = subset(Combined_distribution, Combined_distribution$variable != "Core")

o2 = ggplot(data= Combined_distribution , aes(x=Combined_distribution[,3], y= Combined_distribution[,6], fill = variable)) + geom_point() + xlab("Protein Clusters") + ylab("Pfam Domains") +  geom_smooth(method="lm") + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())+theme(axis.line = element_line(color="black", size = 0.5))

o3 = ggplot(data= Combined_distribution , aes(x=Combined_distribution[,3], y= Combined_distribution[,6], fill = variable)) + geom_point() + xlab("Protein Clusters") + ylab("Pfam Domains") + geom_smooth(method="lm") + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())+theme(axis.line = element_line(color="black", size = 0.5)) + facet_wrap(~ Combined_distribution$variable,scales = c("free"))

multiplot(o2,o3,cols=2)


DF3 = DF1[1:180,]
names(DF3)=c('name','score','value')
DF4 = right_join(DF3, ANVIO_cat,by=('name'))
DF5 = data.frame(family=DF4$group,item=DF4$name,score=DF4$score,value=DF4$value)

polarHistogram(DF5, familyLabel = FALSE)



#----------------------------------------------
#	CORE/Accesory/unique per COG CATEGORY
#----------------------------------------------

COREAnv = AnvioData[,colSums(AnvioData)==60 ]
COREmetadata = ANVIO[ANVIO $protein_cluster_id %in% colnames(COREAnv),]
deduplCOREanvio = distinct(COREmetadata, COG_CATEGORY_ACC, protein_cluster_id)

SINGLAnv = AnvioData[,colSums(AnvioData)==1 ]
SINGLmetadata = ANVIO[ANVIO$protein_cluster_id %in% colnames(SINGLAnv),]
deduplSINGLanvio = distinct(SINGLmetadata, COG_CATEGORY_ACC, protein_cluster_id)

ACCESLAnv = AnvioData[,colSums(AnvioData)>1 & colSums(AnvioData)<60 ]
ACCESmetadata = ANVIO[ANVIO$protein_cluster_id %in% colnames(ACCESLAnv),]
deduplACCESanvio = distinct(ACCESmetadata, COG_CATEGORY_ACC, protein_cluster_id)

ACCESmetadata2 = ACCESmetadata
ACCESmetadata2 $combined = as.factor(paste(ACCESmetadata2 $COG_FUNCTION, ACCESmetadata2 $COG_FUNCTION_ACC,sep='___'))

mat = t(table(droplevels(ACCESmetadata2[c('combined','genome_name')])))

mat2 = t(table(droplevels(ACCESmetadata2[c('COG_FUNCTION_ACC','genome_name')])))
mat2 = as.data.frame.matrix(mat2)

write.csv(mat2,'COG_stamp.txt')

summary(mat2)
cor(mat2)
pca1 = princomp(mat2,scores=TRUE,cor=TRUE) 

pca4 = panpca(mat2, scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16,main="C")

fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 50)
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 100)
fviz_contrib(res.pca, choice = "var", axes = 1:3, top = 100)

fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 30)
fviz_contrib(res.pca, choice = "var", axes = 3, top = 30)
fviz_contrib(res.pca, choice = "var", axes = 4, top = 30)
fviz_contrib(res.pca, choice = "var", axes = 5, top = 30)

plot(pca.res)

STAMP1 = read.table(file = '~/adhaerens_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
data.frame(COG=as.character(STAMP1$COG),p=STAMP1$p.values)
STAMP2 = read.table(file = '~/algicola_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP3 = read.table(file = '~/antarcticus_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP4 = read.table(file = '~/excellens_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP5 = read.table(file = '~/hydrocarb_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP6 = read.table(file = '~/psychro_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP7 = read.table(file = '~/vinifirmus_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP8 = read.table(file = '~/lipolyticus_enriched_list2.tsv', sep = '\t', header = TRUE)

annot_df <- data.frame(adhaerens = subset(STAMP1, COG  %in% COGselect)[,7],algicola = subset(STAMP2, COG  %in% COGselect)[,7],antarcticus = subset(STAMP3, COG  %in% COGselect)[,7])
# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
col = list(adhaerens = circlize::colorRamp2(c(min(annot_df$adhaerens), max(annot_df$adhaerens)),c("lightblue", "purple")),algicola = circlize::colorRamp2(c(min(annot_df$algicola), max(annot_df$algicola)),c("lightblue", "purple")) )
# Create the heatmap annotation
ha <- HeatmapAnnotation(annot_df, col = col)
# Combine the heatmap and the annotation
Heatmap(mat2, name = "mtcars",
        top_annotation = ha)

#NITROGEN RELATED
COGselect_grep = "NITROGEN RELATED"

COGselect = c("COG1116","COG1116-COG4754","COG0715","COG0715-COG2703-COG3706","COG0715-COG3706","COG0600","COG0784-COG0784-COG2198-COG3850-COG4251","COG3062","COG2199-COG2200-COG3850","COG0840-COG3850","COG5013","COG5013-COG5013","COG2180","COG1140","COG3043","COG2181","COG2223","COG2223-COG2223","COG4459","COG3850","COG3850-COG4251","COG3850-COG5001","COG3005","COG3706-COG3850","COG3706-COG3850-COG4564","COG3850","COG1251","COG2146","COG2223-COG2223","COG2223","COG3850-COG4251","COG0243-COG1251","COG3850-COG5001","COG2116","COG2199-COG2200-COG3850","COG1251-COG1773","COG0784-COG0784-COG2198-COG3850-COG4251","COG0457-COG1413-COG3303","COG3706-COG3850","COG0840-COG3850","COG3706-COG3850-COG4564","COG3420","COG3170-COG3420","COG3420-COG3420","COG3420-COG4870","COG4263","COG4314","COG5554","COG2710","COG0840-COG3829-COG5000","COG5000","COG2199-COG2202-COG2203-COG3852-COG5002","COG4191-COG5000-COG5001","COG0347","COG5000-COG5001-COG5002","COG1639-COG2203-COG3852","COG3852-COG5001","COG3706-COG5000","COG2202-COG5000-COG5001","COG3852","COG2203-COG5000-COG5001","COG5000-COG5001","COG0784-COG2205-COG5000","COG0642-COG0745-COG2202-COG2202-COG2203-COG3829-COG5000","COG0784-COG2202-COG2203-COG5000-COG5002","COG0004-COG2202-COG3829-COG5000-COG5001","COG3706-COG3852","COG3197","COG1366-COG3852","COG0745-COG0784-COG2198-COG2202-COG2202-COG3852-COG4251-COG5000","COG0784-COG2202-COG2203-COG2205-COG3614-COG4564-COG5000","COG2199-COG2200-COG3829-COG3852-COG5002","COG0642-COG5000","COG2203-COG3706-COG5000","COG1086-COG2199-COG5000","COG1348","COG0830","COG0378","COG0829","COG2370","COG0804","COG2371","COG0832","COG0831","COG0831-COG0832","COG4413","COG4577","COG4816","COG1881","COG2194","COG4819","COG4812","COG4302","COG4576","COG1881-COG2814","COG3192","COG4810","COG4766","COG4303","COG4820","COG4917","COG4577-COG4577")

#TRANSPORT
COGselect_grep = "TRANSPORT"

COGselect = c("COG1879","COG3842","COG4175","COG0834","COG3127","COG2911-COG5126","COG3090","COG0697","COG4191","COG0531","COG0411-COG4177","COG2059-COG2059","COG1119","COG4664","COG0609","COG1593","COG4521","COG4177","COG3965","COG4536","COG1840","COG0725","COG2911","COG4674","COG1732","COG3181","COG0395","COG2095","COG0614","COG1131","COG0600","COG4166","COG0559","COG3218","COG0474","COG1132","COG0444","COG0823","COG0765","COG1174","COG3113","COG4665","COG1638","COG1120","COG1653","COG3333","COG1116","COG4191-COG5001","COG0811","COG0411","COG3639","COG0577-COG0577","COG4531","COG4143","COG0683","COG4663","COG2067","COG0577","COG4191-COG4191","COG1134","COG0471-COG0471","COG4666","COG1135","COG3833","COG1814","COG0715","COG2391","COG0475-COG1226","COG1117-COG2453","COG1463-COG3008","COG1629","COG0747","COG4191-COG5000-COG5001","COG4590","COG0659","COG4987","COG1117","COG1178","COG2223-COG2223","COG2223","COG2391-COG2391","COG0609-COG0609","COG2202-COG4191","COG3221","COG1173","COG0848","COG1127","COG0003","COG4148","COG0601","COG4149","COG3630","COG2113","COG1172","COG4662","COG4597","COG4114","COG0745-COG2202-COG2205-COG4191","COG1108","COG0517-COG4175","COG4181","COG4107","COG0697-COG3706","COG0811-COG1196","COG0834-COG0834-COG2205","COG1176","COG0842-COG1129","COG1177","COG0842-COG1131-COG1131","COG1126","COG3845","COG0323-COG4191","COG2059","COG4135","COG0842","COG4172","COG2199-COG4191","COG4555","COG0834-COG4191","COG0569","COG3638","COG0784-COG2202-COG2202-COG2203-COG4191","COG0410","COG1637-COG3586","COG1131-COG1131","COG1175","COG0619","COG1726","COG0428","COG2358","COG0569-COG0569","COG0733","COG2984","COG0755","COG0385","COG2216","COG2217-COG2608","COG1118","COG3712","COG2203-COG4191","COG0823-COG1506","COG1463","COG1277","COG2373-COG2911-COG3210-COG4932","COG4607","COG0390","COG0573","COG3839","COG2217","COG1079","COG1121","COG4603","COG0471","COG1123","COG1613","COG2116","COG1292","COG1629-COG4771","COG0697-COG0745-COG4191-COG5002","COG2386","COG2239-COG4175","COG0003-COG0003","COG2217-COG2608-COG2608","COG0784-COG3437-COG4191","COG0410-COG3945","COG2197-COG2202-COG2202-COG4191-COG4936","COG1125","COG2871","COG4558","COG0834-COG0834-COG0834-COG3706","COG1742","COG4150","COG4598","COG0834-COG0834-COG3706","COG0517-COG1125","COG0697-COG2510","COG4619","COG0471-COG0471-COG0569","COG0310","COG0488","COG2998","COG4176","COG0598","COG4535","COG0697-COG1305","COG3639-COG3639","COG3201","COG4591","COG0834-COG0834","COG1629-COG4772","COG4160","COG3069","COG2060","COG0823-COG0823-COG0841-COG1361-COG1572-COG1572-COG1657-COG4733","COG2203-COG4191-COG5001","COG0471-COG0471-COG0569-COG3273","COG4988","COG3470","COG4592","COG0226","COG3703","COG0784-COG2202-COG2202-COG4191","COG3840","COG1805-COG1805","COG0168","COG0591-COG4191","COG2270","COG1129","COG1129-COG4177","COG1288","COG2869","COG2854","COG4133","COG5295","COG2076","COG1347","COG2911-COG5295","COG0471-COG0471-COG0490","COG2156","COG0767","COG0784-COG4191","COG4608","COG4134","COG0823-COG0823-COG1506","COG0318-COG0697","COG4239","COG0842-COG1511","COG4215","COG1883","COG0697-COG0784-COG2198-COG2199-COG2200-COG2205","COG0834-COG2205-COG3437","COG1277-COG3225-COG3225","COG2056","COG3706-COG4191","COG3192","COG1116-COG4754","COG0715-COG2703-COG3706","COG4191-COG4192","COG0053","COG1918","COG1479-COG3586","COG0842-COG1131","COG2199-COG3221-COG5002","COG0577-COG0577-COG1511","COG1268","COG4559","COG2209","COG0515-COG2203-COG3899-COG4191","COG1226-COG4651","COG0581","COG0834-COG3706","COG2853","COG0370","COG0025-COG0569","COG4572","COG1122","COG2211","COG0715-COG3706","COG0747-COG1305","COG0515-COG0683","COG2011","COG3322-COG4192-COG4251","COG2911-COG3210-COG3210","COG0475-COG0490-COG1226","COG3437-COG3829-COG4191-COG5001","COG1464","COG0555","COG0823-COG3227","COG0475","COG4618","COG3437-COG4191","COG1687","COG0226-COG2885","COG0683-COG3290-COG5001","COG2217-COG3350","COG2931-COG2931-COG5295","COG0577-COG1136","COG0444-COG4608","COG4779","COG3158","COG2239","COG2031","COG1914","COG2911-COG2982","COG4594","COG4174","COG0559-COG1413","COG3558","COG4176-COG4176","COG4137","COG2911-COG3210-COG3210-COG3210-COG3595","COG3447-COG4191-COG5001","COG1101","COG2202-COG3221-COG5001-COG5002","COG4120","COG1593-COG3090","COG1120-COG3593","COG4413","COG2967","COG2199-COG2703-COG4191","COG5265","COG0823-COG5492-COG5492-COG5492","COG1629-COG4774","COG4208","COG0784-COG2202-COG4191","COG2202-COG2202-COG4191-COG4251","COG1119-COG1119","COG4136")

#phosphonate
COGselect_grep = "phosphonate"

COGselect = c("COG3626","COG3639","COG3221","COG4107","COG3638","COG3627","COG3454","COG3639-COG3639","COG3625","COG2199-COG3221-COG5002","COG3624","COG4778","COG2202-COG3221-COG5001-COG5002")

#cytochrome
COGselect_grep = "cytochrome"

COGselect = c("COG2124","COG3005","COG3125","COG3278","COG2193","COG3245","COG4987","COG2010-COG4993","COG1858","COG2009","COG2010","COG1622-COG2132","COG2010-COG2010","COG3336","COG4736","COG0843","COG3258","COG2863","COG2863-COG2863-COG3258","COG1622-COG2010","COG3474","COG3909","COG2133-COG3291-COG3828-COG4654","COG1138","COG2864","COG4235","COG3088","COG0755","COG1845","COG2993-COG3278","COG2010-COG4733","COG4117","COG2863-COG2863","COG2386","COG1622","COG3038","COG1999","COG0785","COG2332","COG3658","COG4988","COG4654","COG1413-COG2010-COG2133","COG4133","COG2993","COG2010-COG2010-COG2010","COG3175","COG0369-COG2124","COG0457-COG1413-COG3303","COG2010-COG3391-COG3391","COG1858-COG3391-COG3391","COG1858-COG1858-COG1858","COG3043","COG2010-COG4319","COG1294","COG5274","COG1271","COG1018-COG2124","COG3346","COG0843-COG1845","COG2857","COG1858-COG3291-COG3391-COG3391","COG1290")


#ATP synthase
COGselect_grep = "ATP synthase"

COGselect = c("COG0224","COG0355","COG0055","COG0356","COG0711-COG0712","COG3312","COG0636","COG0711","COG0056","COG0712")

#sulfur associated
COGselect_grep = "sulfur"


COGselect = c("COG2168","COG1218","COG0600","COG2141","COG3119","COG0607","COG0225","COG1651","COG4232","COG1416","COG1116","COG0175","COG1104","COG0526","COG0715","COG2391","COG0659","COG1054","COG2391-COG2391","COG0175-COG3415","COG1553","COG1262","COG0446","COG2897","COG1331-COG4233","COG3011","COG2151","COG0229","COG2923","COG2836","COG0369","COG0491-COG0607","COG4424","COG0425","COG2151-COG3677","COG1262-COG4976","COG0369-COG3182","COG1118","COG0520","COG2761","COG0404-COG0446","COG0641","COG4117","COG1613","COG0526-COG0682","COG0155","COG2846","COG0446-COG1902","COG1262-COG1409","COG4150","COG1145-COG1148","COG2015","COG0491-COG0607-COG0607","COG1495","COG0225-COG0229","COG2895","COG0529-COG2895","COG0489-COG2151","COG4312","COG0476-COG0607-COG1977","COG0306","COG2104","COG0306-COG0306","COG0446-COG1773","COG1116-COG4754","COG0715-COG2703-COG3706","COG2046","COG0369-COG2124","COG1262-COG5635","COG0607-COG0640","COG4630","COG0446-COG3453","COG0301-COG0607","COG0715-COG3706","COG0555","COG2920","COG3531","COG2873","COG0607-COG0607-COG0607-COG0607","COG0607-COG0607-COG2897","COG0607-COG2897-COG2897","COG2166","COG1262-COG1352","COG1262-COG2227","COG4208","COG0529")

#cobalamin
COGselect_grep = "cobalamin"

COGselect = c("COG4206","COG1120","COG2885-COG4547-COG5384-COG5644","COG4547","COG0620","COG4206-COG4771","COG4812","COG2185","COG0620-COG0620","COG1270","COG0368","COG1429","COG0646-COG1410","COG0646","COG1410","COG1429-COG1429","COG1120-COG3593","COG1703-COG1884-COG2185")

#palsmid
COGselect_grep = "plasmid"

COGselect = c("COG4643","COG5527","COG4691","COG3093","COG3378","COG2856-COG3093","COG4643-COG5519","COG5302","COG5655","COG3668","COG3030","COG4643-COG4643-COG5519","COG3549","COG3657","COG5534","COG4227-COG4643","COG0358-COG3378")

#transpos
COGselect_grep = "transpos"

COGselect = c("COG2963","COG1943","COG3666","COG4584","COG3039","COG3677","COG5421","COG2801","COG3328","COG3464","COG3385","COG3415","COG3436","COG0175-COG3415","COG3293","COG3547","COG2826","COG5659","COG4644","COG3385-COG5659","COG2151-COG3677","COG5464","COG3316","COG5433","COG0675","COG3335-COG3415","COG3335","COG2842","COG2801-COG2963-COG2963","COG2801-COG2963","COG2801-COG3415","COG4565-COG4584","COG2771-COG4584","COG3666-COG3666","COG1398-COG3464","COG3385-COG5421","COG2963-COG2963","COG3293-COG5421","COG1343","COG1857","COG1203","COG1468","COG3513","COG1518","COG3649")

#CRIPSR
COGselect_grep = "CRIPSR"

COGselect = c("COG1343","COG1857","COG1203","COG1468","COG3513","COG1518","COG3649")

#integrase
COGselect_grep = "integrase"

COGselect = c("COG0582","COG0582-COG4974")

#nitrous

COGselect = c("COG3420","COG3170-COG3420","COG3420-COG3420","COG3420-COG4870","COG4263","COG4314","COG0742","COG0348-COG3901","COG3256","COG4548")

#pilus
COGselect_grep = "pilus"

COGselect = c("COG4970","COG4967","COG4962","COG3419","COG4961","COG4726","COG2804","COG4969","COG5010","COG4968","COG4963","COG4966","COG2207-COG3063-COG5616","COG3170","COG3170-COG3420","COG1530-COG3170","COG3745","COG3167","COG3166","COG3063","COG3170-COG4313","COG0501-COG4969","COG3121","COG4964","COG0457-COG3063","COG3166-COG4972","COG0249-COG0553-COG0827-COG3170-COG4646","COG0643-COG2198-COG3170","COG2304-COG3419","COG3168","COG3215","COG2304-COG2304-COG3419","COG5008","COG4965","COG4972","COG2518-COG2804","COG2064","COG2805","COG3706-COG4963","COG3152-COG4969","COG0457-COG3170","COG2804-COG3437","COG3847")
#secretiuon
COGselect = c("COG2831","COG4104","COG3267","COG2804","COG3520","COG3157","COG2304","COG3266-COG3267","COG3211","COG4959","COG3843","COG2257","COG3501","COG3701","COG1886","COG3451","COG1459","COG2165","COG3149","COG4796","COG3505","COG2304-COG2931","COG3517","COG5621","COG3519","COG5501","COG4964","COG1357-COG1357-COG2948","COG1450","COG5268","COG3838","COG3846","COG1157","COG3522","COG2948","COG0805","COG0457-COG2304","COG3156","COG1989","COG1025","COG2304-COG3419","COG1766","COG1711-COG2304","COG2304-COG2304-COG3419","COG3031","COG3523","COG2340-COG5640","COG3504","COG3521","COG0457-COG4796","COG2503-COG3568","COG1462","COG1729-COG2304","COG5510","COG2518-COG2804","COG3516","COG2304-COG2885","COG1657-COG2304","COG3515-COG5283-COG5412","COG1317","COG0084","COG3455","COG2804-COG3437","COG3515","COG3418","COG1450-COG4796","COG3518","COG5497","COG1716-COG3456","COG4795","COG3211-COG3386","COG1409-COG2831","COG3297")

#phage
COGselect_grep = "phage"

COGselect = c("COG5511","COG5525","COG4220","COG4643","COG3497","COG1842","COG3645","COG4540","COG4387","COG5004","COG2369","COG4383","COG3378","COG4653","COG3636","COG2932","COG3501","COG5005","COG4228","COG4397","COG4643-COG5519","COG5301","COG5283","COG4733","COG3778","COG4679","COG3772","COG4518","COG4388","COG1372-COG2369","COG3948","COG5280","COG1361-COG4733","COG4695","COG3617","COG4382","COG5471","COG4197","COG4373","COG2915","COG2369-COG4695","COG4385","COG3299","COG3654-COG3943","COG3498","COG2010-COG4733","COG3499","COG3561","COG4626","COG4678","COG4733-COG4733","COG5003","COG4823","COG2842","COG3728","COG1396-COG2932","COG1511","COG4384","COG3210-COG4733-COG4733-COG4733-COG4733-COG4733-COG4733-COG4733-COG4733-COG4733-COG4733-COG4733-COG4733-COG4733-COG4733","COG4733-COG4934","COG5377","COG3030","COG4381","COG0419-COG3941","COG3409-COG3772","COG0741-COG3941","COG0823-COG0823-COG0841-COG1361-COG1572-COG1572-COG1657-COG4733","COG3941","COG3500","COG4643-COG4643-COG5519","COG4379","COG1196-COG5283","COG5281","COG4416","COG4580","COG1061-COG2227-COG2932-COG3886","COG3654","COG3646","COG3628","COG1783","COG3747","COG0842-COG1511","COG4675-COG5301","COG3740","COG5362","COG4718","COG0577-COG0577-COG1511","COG2931-COG2931-COG4733","COG4422","COG5362-COG5410","COG4570","COG3515-COG5283-COG5412","COG3600","COG5614","COG3291-COG3291-COG4733-COG4733","COG3497-COG3497","COG5518","COG1983","COG1842-COG4656","COG4672","COG4386","COG4227-COG4643","COG0358-COG3378")

#AROMATIC
COGselect_grep = "AROMATIC"

COGselect = c("COG3384","COG0079","COG4664","COG4665","COG4663","COG1448","COG3805","COG2368")

#Chemotaxis
COGselect_grep = "Chemotaxis"

COGselect = c("COG0840","COG0784-COG2198-COG2202-COG2205-COG3290-COG3706-COG3706-COG3829","COG0840-COG2202","COG2201","COG0784-COG0835","COG0835","COG0840-COG3829","COG0840-COG3829-COG5000","COG0784-COG2198-COG3614-COG5002","COG0784-COG0784-COG2198-COG2205","COG0784","COG0784-COG0784-COG2198-COG2205-COG3437-COG5001","COG0643-COG0784-COG2198-COG2198-COG2198-COG2198-COG2198","COG0643-COG0745-COG2198-COG2198-COG2198-COG2198-COG2198-COG2198","COG0643-COG0784-COG2198-COG2198-COG2198-COG2198-COG2198-COG2198","COG0745-COG0784-COG2198-COG5002","COG0784-COG2205-COG2703","COG0784-COG2205-COG4251","COG0784-COG2198-COG5002","COG0840-COG3290","COG0784-COG4251","COG0784-COG2202-COG2202-COG2203-COG4191","COG0642-COG0784-COG2198-COG2770","COG1871","COG0784-COG0784-COG2198-COG2202-COG2205-COG2461","COG0840-COG4564","COG0642-COG0784-COG0784-COG2198-COG3447-COG3614","COG1196-COG1352-COG2201-COG3829-COG3829-COG5001","COG1352-COG2201-COG3829-COG4717-COG5001","COG0643-COG2198","COG0642-COG0784-COG2202","COG0784-COG3437-COG4191","COG0840-COG5278","COG0840-COG5278-COG5278","COG0775-COG0784","COG0784-COG0784-COG2198-COG3452-COG5002","COG0784-COG2205-COG5000","COG0643-COG2198-COG3170","COG0784-COG2202-COG2203-COG5000-COG5002","COG0784-COG2203-COG2203-COG5001","COG0840-COG2202-COG2202","COG1352","COG0784-COG2200-COG4566","COG1776-COG2201","COG0784-COG2202-COG2202-COG4191","COG0784-COG2205","COG0784-COG2202-COG2203-COG2205","COG0835-COG3706","COG0784-COG5002","COG0784-COG4191","COG0784-COG2202-COG2205","COG0457-COG0784","COG0784-COG0784-COG2198-COG5002","COG0784-COG2198-COG2202-COG5001","COG0784-COG0784-COG2198-COG3850-COG4251","COG0697-COG0784-COG2198-COG2199-COG2200-COG2205","COG0745-COG0784-COG2198-COG2202-COG2202-COG3852-COG4251-COG5000","COG0784-COG0784-COG2198-COG2202-COG2202-COG2205-COG3829-COG5002","COG0784-COG2202-COG2203-COG2205-COG3614-COG4564-COG5000","COG0784-COG2202-COG2202-COG3706","COG0784-COG2197-COG2205","COG0840-COG2202-COG2202-COG2202","COG0784-COG2198-COG2202-COG2205-COG3614","COG0835-COG3266","COG0745-COG0784-COG2198-COG4251","COG0784-COG2197-COG2198-COG2205","COG0784-COG3290-COG3706","COG3143","COG0784-COG2199-COG2200","COG0591-COG0784-COG2205","COG0840-COG1893","COG0840-COG4372","COG0784-COG2205-COG3300","COG1262-COG1352","COG0784-COG2198-COG2205","COG0840-COG3850","COG0643","COG0643-COG0643-COG0745-COG2198","COG1196-COG1352-COG2201","COG0784-COG2205-COG3726","COG0784-COG2202-COG4191","COG0784-COG2202-COG2202-COG2202-COG2205","COG1352-COG2201-COG3352","COG1352-COG2201-COG2461")

#Iron, Ferri, Fe
COGselect_grep = "Iron, Ferri, Fe"

COGselect = c("COG0609","COG1120","COG4264","COG0609-COG0609","COG4615","COG2977","COG2375","COG1120-COG3593","COG0007-COG1648","COG2837","COG0276","COG1648","COG2138","COG1629","COG1629-COG4771","COG1629-COG4772","COG1629-COG4774","COG3182","COG2710","COG0316","COG4114","COG0369-COG3182","COG3712","COG0390","COG3016","COG2846","COG4619","COG3487","COG5007","COG3022","COG4943","COG0316-COG0694","COG4630","COG0371")

#carbon
COGselect = c("COG0600","COG0663","COG1966","COG1116","COG0715","COG1551","COG3427","COG2080","COG1529-COG3427","COG1719-COG3829","COG0288","COG1116-COG4754","COG0715-COG2703-COG3706","COG0715-COG3706","COG1719","COG2080-COG3427")

#carbonic anhydrase and cyanate lyase
COGselect = c("COG0288","COG2807","COG1513")

#glyco
COGselect = c("COG0463","COG0438","COG1216","COG1819","COG1215","COG0741","COG0438-COG0457-COG0457-COG1216","COG0366","COG0438-COG3072","COG0438-COG0463","COG2951","COG1807","COG1472","COG0297","COG3173","COG1331-COG4233","COG1216-COG1216","COG1397","COG0546","COG3306","COG0122","COG0266","COG0438-COG1409","COG0438-COG1216-COG2227","COG1573","COG1194","COG0463-COG3307","COG3194","COG3980","COG0438-COG0613","COG0463-COG1215","COG0030","COG0463-COG0500-COG1215","COG1216-COG1216-COG2227","COG0692","COG0438-COG0707","COG0438-COG0561","COG0645-COG2187","COG0438-COG0438-COG0438","COG3408","COG4671","COG0741-COG0791","COG0438-COG0457-COG0457-COG0457-COG0457-COG0457-COG0457-COG0859-COG0859-COG0859-COG0859","COG1397-COG1397-COG1397","COG1817","COG5520","COG2818","COG1442","COG0741-COG3941","COG0438-COG0457-COG1216-COG1216","COG0296-COG0366","COG2055","COG1215-COG1216","COG0615-COG1216","COG1523","COG3222","COG0438-COG0463-COG1216-COG2227","COG0438-COG1216-COG1216","COG0122-COG2169","COG0500-COG1216","COG2261","COG1216-COG1887","COG1059","COG0438-COG2226-COG3587","COG1215-COG1215","COG2951-COG3409","COG0366-COG3280","COG0438-COG1216-COG2226-COG2835","COG4623","COG2943","COG0741-COG1388-COG1388-COG1388","COG3178","COG0438-COG0535","COG2187","COG1215-COG2246","COG0438-COG0438","COG2152","COG1501","COG1216-COG4261","COG1215-COG4261","COG2273-COG3306","COG0546-COG0546")


#tRNA
COGselect = c("COG0324","COG0590","COG0215","COG0533","COG1639-COG2606","COG0154","COG0252","COG3642","COG2265","COG1720","COG0219","COG0073","COG1690","COG2872","COG0223-COG1670","COG0162","COG1234","COG0173","COG2603","COG0124","COG0442","COG0486","COG2606","COG1190","COG0008","COG1444","COG0060","COG0064","COG0073-COG0143","COG0172","COG0223","COG0566","COG0180","COG0013","COG0037","COG0018","COG0564","COG1921","COG0172-COG5653","COG0751","COG0354","COG0343","COG0665-COG4121","COG0193","COG4123","COG0441","COG0820","COG0525","COG1179","COG0445","COG4445","COG0220","COG4694","COG0017","COG0009","COG1214","COG1746","COG0072-COG0073","COG0617","COG0373","COG0042","COG0223-COG0726","COG0495-COG0495","COG0101","COG2935","COG0621","COG0585","COG0301-COG0607","COG0565","COG0495","COG0752","COG1490","COG0482","COG0809","COG0336","COG0016","COG0802","COG2520","COG2360","COG0130","COG1190-COG2512","COG0721")

#aerobic
COGselect = c("COG2080","COG0243-COG1251","COG0243-COG2906","COG1328","COG0243","COG0243-COG3383","COG0243-COG0437","COG2080-COG3427")


#nitrate_reductase
COGselect_grep = "nitrate_reductase"

COGselect = c("COG0324","COG3062","COG1140","COG2180","COG5013","COG5013-COG5013","COG2181","COG3043","COG4459")


#anaplerotic
COGselect_grep = "anaplerotic"

COGselect = c("COG5524","COG0288","COG1866","COG0458","COG0505","COG0600","COG1116","COG0715","COG2352","COG0439","COG0511-COG5016","COG0511")

#quinone
COGselect_grep = "quinone"

COGselect = c("COG1034","COG0713","COG2226","COG0604","COG2941","COG1007","COG0661","COG0377","COG1009","COG5424","COG1008","COG0435","COG0457-COG0457-COG2226-COG4976","COG1726","COG0789-COG2226","COG0838","COG0839","COG2249","COG0236-COG0300-COG0604-COG3321","COG2871","COG1005","COG0649-COG0852","COG3165","COG1009-COG2111","COG1805-COG1805","COG2869","COG1347","COG3865","COG0438-COG2226-COG3587","COG1143","COG2209","COG0438-COG1216-COG2226-COG2835","COG1905","COG0674-COG1013-COG1014-COG1143","COG1009-COG2111-COG2111","COG1894-COG1905","COG0382-COG2249","COG1894")


#repair
COGselect_grep = "repair"

COGselect = c("COG0419-COG3593","COG0420","COG1195-COG3593","COG1119","COG0419","COG1273","COG1690","COG0323","COG3145","COG1197","COG0389","COG0497-COG1061","COG0323-COG4191","COG2189-COG2852","COG3723","COG0497","COG0536","COG1533","COG0354","COG3727","COG2846","COG0249-COG0553-COG0827-COG3170-COG4646","COG0419-COG3941","COG0497-COG3593","COG2852","COG1061-COG1112-COG2852","COG0610-COG0610-COG2852","COG2924","COG0353","COG0286-COG0497","COG2852-COG5340","COG1195","COG0323-COG2205","COG0249","COG0062-COG0063","COG1112-COG1199-COG2852","COG1381","COG0419-COG2274","COG0419-COG3950","COG2003","COG0323-COG1413","COG0507-COG1112-COG2852","COG2354","COG3695","COG1119-COG1119","COG0419-COG1195")


#malate synthase
COGselect = c("COG2225","COG2224")

#efflux
COGselect_grep = "efflux"

COGselect = c("COG0841","COG2814","COG3965","COG3696","COG1566","COG0845","COG1280","COG1230","COG0845-COG5569","COG0534","COG0845-COG1994","COG0823-COG0823-COG0841-COG1361-COG1572-COG1572-COG1657-COG4733","COG5569","COG0845-COG3667-COG5569","COG1230-COG2608","COG3696-COG3696","COG1668","COG1881-COG2814","COG0798","COG1971","COG0803-COG0845","COG1322-COG1566","COG0204-COG2814","COG1538-COG3696")


#metal efflux
COGselect_grep = "metalefflux"

COGselect = c("COG3965","COG3696","COG1230","COG1230-COG2608","COG3696-COG3696","COG1538-COG3696","COG1814","COG1108","COG1121","COG1971","COG1914")

#Beta-lactamase
COGselect_grep = "Beta-lactamase"

COGselect = c("COG0491","COG2333","COG2367","COG0658-COG2333","COG1236","COG0491-COG0607","COG2220","COG1680","COG1680-COG3653","COG2015","COG0491-COG0607-COG0607","COG2248","COG2602-COG3904","COG3725","COG0810-COG4219","COG1396-COG1680-COG3296","COG0491-COG0494")

#UREA
COGselect_grep = "UREA"

COGselect = c("COG1984","COG0439-COG0511-COG1984-COG2049","COG0830","COG0378","COG0829","COG2370","COG0804","COG2371","COG0832","COG0831","COG0831-COG0832","COG4413")

#glycosid
COGselect = c("COG0366","COG1472","COG3173","COG1331-COG4233","COG0645-COG2187","COG0296-COG0366","COG0366-COG3280","COG3178","COG2187")

#glyco
COGselect_grep = "glyco"

COGselect = c("COG0463","COG0438","COG1216","COG1819","COG1215","COG0741","COG0438-COG0457-COG0457-COG1216","COG0366","COG0438-COG3072","COG0438-COG0463","COG2951","COG1807","COG1472","COG0297","COG3173","COG1331-COG4233","COG1216-COG1216","COG1397","COG0546","COG3306","COG0122","COG0266","COG0438-COG1409","COG0438-COG1216-COG2227","COG1573","COG1194","COG0463-COG3307","COG3194","COG3980","COG0438-COG0613","COG0463-COG1215","COG0030","COG0463-COG0500-COG1215","COG1216-COG1216-COG2227","COG0692","COG0438-COG0707","COG0438-COG0561","COG0645-COG2187","COG0438-COG0438-COG0438","COG3408","COG4671","COG0741-COG0791","COG0438-COG0457-COG0457-COG0457-COG0457-COG0457-COG0457-COG0859-COG0859-COG0859-COG0859","COG1397-COG1397-COG1397","COG1817","COG5520","COG2818","COG1442","COG0741-COG3941","COG0438-COG0457-COG1216-COG1216","COG0296-COG0366","COG2055","COG1215-COG1216","COG0615-COG1216","COG1523","COG3222","COG0438-COG0463-COG1216-COG2227","COG0438-COG1216-COG1216","COG0122-COG2169","COG0500-COG1216","COG2261","COG1216-COG1887","COG1059","COG0438-COG2226-COG3587","COG1215-COG1215","COG2951-COG3409","COG0366-COG3280","COG0438-COG1216-COG2226-COG2835","COG4623","COG2943","COG0741-COG1388-COG1388-COG1388","COG3178","COG0438-COG0535","COG2187","COG1215-COG2246","COG0438-COG0438","COG2152","COG1501","COG1216-COG4261","COG1215-COG4261","COG2273-COG3306","COG0546-COG0546")

#fructose
COGselect_grep = "fructose"

COGselect = c("COG0483","COG1080-COG1925-COG4668","COG3588","COG1299-COG1445-COG1445","COG1830","COG1105","COG0191","COG1762")

#mannitol
COGselect_grep = "mannitol"

COGselect = c("COG0483","COG1080-COG1925-COG4668","COG3588","COG1299-COG1445-COG1445","COG1830","COG1105","COG0191","COG1762","COG0246","COG4664","COG4665","COG4663","COG1080-COG1925-COG4668","COG1762")
#acetate
COGselect_grep = "acetate"

COGselect = c('COG1788','COG1541','COG0282','COG1060','COG2057','COG3630','COG0511-COG5016','COG4670','COG4147','COG1883','COG4689','COG3970')

#light
COGselect_grep = "light"

COGselect = c('COG0511-COG5016','COG3829-COG4251','COG2203-COG3300-COG4251','COG3850-COG4251','COG2202-COG4251-COG5001','COG0784-COG2205-COG4251','COG0784-COG4251','COG3614-COG4251','COG3706-COG4251','COG2202-COG2202-COG2202-COG3829-COG4251','COG0745-COG2197-COG2202-COG2202-COG3829-COG4251','COG0784-COG0784-COG2198-COG3850-COG4251','COG0745-COG0784-COG2198-COG2202-COG2202-COG3852-COG4251-COG5000','COG2202-COG2202-COG2202-COG2202-COG2203-COG2203-COG3290-COG3290-COG3290-COG4251','COG4251','COG0745-COG0784-COG2198-COG4251','COG2202-COG4251','COG2202-COG3829-COG4251-COG5002','COG3322-COG4192-COG4251','COG2202-COG2202-COG4191-COG4251')

#repair
COGselect_grep = "repair"

COGselect = c('COG0419-COG3593','COG0420','COG1195-COG3593','COG1119','COG0419','COG1273','COG1690','COG0323','COG3145','COG1197','COG0389','COG0497-COG1061','COG0323-COG4191','COG2189-COG2852','COG3723','COG0497','COG0536','COG1533','COG0354','COG3727','COG2846','COG0249-COG0553-COG0827-COG3170-COG4646','COG0419-COG3941','COG0497-COG3593','COG2852','COG1061-COG1112-COG2852','COG0610-COG0610-COG2852','COG2924','COG0353','COG0286-COG0497','COG2852-COG5340','COG1195','COG0323-COG2205','COG0249','COG0062-COG0063','COG1112-COG1199-COG2852','COG1381','COG0419-COG2274','COG0419-COG3950','COG2003','COG0323-COG1413','COG0507-COG1112-COG2852','COG2354','COG3695','COG1119-COG1119','COG0419-COG1195')
#photo
COGselect = c('COG1119','COG4447','COG0415','COG3046','COG1533','COG4447-COG4447','COG4447-COG4447-COG4447','COG1119-COG1119','COG0783')


#stress
COGselect_grep = "stress"

COGselect = c('COG0589','COG2334','COG2310','COG0589-COG0589','COG2310-COG4110','COG3861','COG3861-COG3861','COG3871','COG0271','COG5007','COG2976','COG1825','COG0783','COG1983','COG1217')
#superoxide dismutase and catalase
COGselect_grep = "superoxide dismutase and catalase"

COGselect = c('COG0753','COG0693-COG0753','COG0376','COG0605','COG2032')
#peroxiredoxin, thioredoxin
COGselect_grep = "peroxiredoxin, thioredoxin"

COGselect = c('COG0526','COG1225','COG0450','COG2044','COG0678','COG3118','COG0526-COG0682','COG0492','COG2077','COG0678-COG0695','COG2143')
#redoxin
COGselect = c('COG0543-COG0633','COG0348','COG1018','COG0526','COG1225','COG4454','COG0348-COG3901','COG4656-COG4656','COG0450','COG2146','COG2132','COG2132-COG2132','COG0695','COG0633','COG1014-COG4231','COG2044','COG2878','COG1331-COG4233','COG1145','COG0543-COG1018','COG1622-COG2132','COG0426','COG0678','COG4659','COG3118','COG0243-COG2906','COG0543-COG1018-COG3239','COG2906','COG0526-COG0682','COG1018-COG3576-COG3576','COG0492','COG1018-COG1018','COG1018-COG3000','COG1145-COG1149','COG1145-COG1148','COG1139','COG1017-COG1018','COG0247-COG0277-COG1145','COG4660','COG2077','COG4658','COG4657','COG1393','COG1773','COG1773-COG1773','COG0446-COG1773','COG1251-COG1773','COG4656','COG0493-COG1144','COG0278','COG0678-COG0695','COG0674-COG1013-COG1014-COG1143','COG0674-COG1013-COG1014-COG1146','COG2143','COG0644-COG2440','COG1018-COG2124','COG1842-COG4656','COG1196-COG4656','COG1141','COG1018-COG3576')

#superoxide dismutase and catalase, peroxiredoxin, thioredoxin
COGselect_grep = "superoxide dismutase and catalase, peroxiredoxin, thioredoxin"

COGselect = c('COG0753','COG0693-COG0753','COG0376','COG0605','COG2032','COG0526','COG1225','COG0450','COG2044','COG0678','COG3118','COG0526-COG0682','COG0492','COG2077','COG0678-COG0695','COG2143')

#amino acid permease
COGselect = c('COG0411-COG4177','COG4177','COG0559','COG0765','COG1296','COG0814','COG4597','COG1129-COG4177','COG0833','COG0559-COG1413')
#amino acid ABC
COGselect_grep = "amino-acid transport"

COGselect = c('COG0834','COG0411-COG4177','COG4177','COG0559','COG0765','COG0411','COG0683','COG4597','COG0834-COG0834-COG2205','COG1126','COG0834-COG4191','COG0410','COG0410-COG3945','COG0834-COG0834-COG0834-COG3706','COG0834-COG0834-COG3706','COG0834-COG0834','COG1129-COG4177','COG0834-COG2205-COG3437','COG0834-COG3706','COG0515-COG0683','COG0683-COG3290-COG5001','COG0559-COG1413')

#grep amino_acid
COGselect_grep = "amino-acid"

COGselect = c('COG0834','COG0531','COG0411-COG4177','COG4177','COG2095','COG0559','COG0765','COG1296','COG0411','COG0683','COG0814','COG0115','COG0665','COG4597','COG0834-COG0834-COG2205','COG1126','COG0834-COG4191','COG0410','COG1247','COG1042-COG1247','COG0410-COG3945','COG0665-COG4121','COG0834-COG0834-COG0834-COG3706','COG5322','COG0388-COG1247','COG0834-COG0834-COG3706','COG0834-COG0834','COG0404-COG0665','COG1129-COG4177','COG0834-COG2205-COG3437','COG0833','COG0834-COG3706','COG0515-COG0683','COG1687','COG0683-COG3290-COG5001','COG0560-COG3830','COG0559-COG1413')

#grep branched amino_acid
COGselect_grep = "branched amino_acid"

COGselect = c('COG0411-COG4177','COG4177','COG0559','COG1296','COG0411','COG0683','COG0115','COG0410','COG0410-COG3945','COG1129-COG4177','COG0515-COG0683','COG1687','COG0683-COG3290-COG5001','COG0559-COG1413')

#grep metal resitance
COGselect_grep = "metal resitance"

COGselect = c('COG1055','COG1393','COG0798','COG0053','COG3667','COG3015','COG3015-COG3187','COG0845-COG3667-COG5569','COG2132','COG2132-COG2132','COG3793','COG4103','COG1275','COG2310-COG4110','COG0861','COG0500-COG3615','COG0484-COG4103','COG3615','COG3853','COG0428','COG5446','COG3965','COG1230','COG1230-COG2608','COG3696','COG3696-COG3696','COG1538-COG3696','COG1814','COG1108','COG1121','COG1971','COG1914')

#arginine
COGselect_grep = "arginine"

COGselect = c('COG3138','COG1586','COG2235','COG4160','COG1279','COG4215','COG2957','COG1166','COG3724','COG3880')
#GGDEF
COGselect_grep = "GGDEF"

COGselect = c('COG0784-COG2198-COG2202-COG2205-COG3290-COG3706-COG3706-COG3829','COG2199-COG2200','COG3437-COG5001','COG0745-COG2198-COG3706','COG3437-COG3829-COG5001','COG3706-COG5002','COG5001','COG3706','COG0745-COG3290-COG5001','COG2200-COG3706','COG2202-COG2202-COG3706','COG4191-COG5001','COG4564-COG5001','COG0517-COG2199-COG2200','COG2199-COG2203','COG2202-COG2202-COG5001','COG3706-COG5001','COG2203-COG5001','COG2199','COG2199-COG2202-COG2203-COG3852-COG5002','COG0745-COG5001','COG0784-COG0784-COG2198-COG2205-COG3437-COG5001','COG2202-COG5001','COG4191-COG5000-COG5001','COG2905-COG2905-COG5001','COG2199-COG2200-COG3287','COG2199-COG2203-COG3829','COG3706-COG3829','COG2202-COG3614-COG5001','COG5000-COG5001-COG5002','COG2199-COG2703','COG0697-COG3706','COG2199-COG2200-COG2204','COG3322-COG3706','COG3852-COG5001','COG2199-COG4191','COG3706-COG5000','COG2202-COG4251-COG5001','COG2202-COG3706','COG2202-COG5000-COG5001','COG3447-COG3829-COG5001','COG2203-COG3300-COG3706','COG2199-COG2905-COG2905','COG2202-COG3706-COG5001','COG5001-COG5002','COG3850-COG5001','COG2202-COG3284-COG5001','COG2199-COG2200-COG2203','COG3829-COG5001','COG3829-COG5001-COG5002','COG1196-COG1352-COG2201-COG3829-COG3829-COG5001','COG0745-COG2202-COG5001-COG5002','COG1352-COG2201-COG3829-COG4717-COG5001','COG2202-COG2202-COG2203-COG5001','COG2203-COG5000-COG5001','COG2203-COG3706','COG3452-COG5001','COG2203-COG2203-COG5001','COG3706-COG4251','COG2203-COG2203-COG2203-COG2203-COG2203-COG5001','COG1716-COG5001','COG2205-COG3829-COG5001','COG5000-COG5001','COG2199-COG2200-COG2203-COG3290','COG0834-COG0834-COG0834-COG3706','COG2202-COG2203-COG3290-COG5001','COG2202-COG2203-COG5001','COG0834-COG0834-COG3706','COG2199-COG2200-COG3279','COG0457-COG2202-COG5001','COG2202-COG2461-COG5001','COG0004-COG2202-COG3829-COG5000-COG5001','COG3284-COG5001','COG0784-COG2203-COG2203-COG5001','COG2199-COG3300','COG2203-COG4191-COG5001','COG2199-COG2200-COG3290','COG3706-COG3852','COG2200-COG2204-COG3706','COG2203-COG3290-COG5001','COG1835-COG2199','COG2199-COG3829','COG0457-COG0515-COG2203-COG2203-COG3706-COG3899','COG0835-COG3706','COG2198-COG3706','COG2202-COG3279-COG5001','COG0784-COG2198-COG2202-COG5001','COG0745-COG3706','COG2199-COG2200-COG3850','COG2199-COG3447','COG3264-COG3706','COG0697-COG0784-COG2198-COG2199-COG2200-COG2205','COG3706-COG4191','COG0715-COG2703-COG3706','COG0784-COG2202-COG2202-COG3706','COG2208-COG3706','COG3447-COG3706','COG2199-COG3221-COG5002','COG2203-COG3829-COG5001','COG2202-COG2202-COG2203-COG2203-COG2203-COG5001','COG0834-COG3706','COG3614-COG5001','COG3706-COG4963','COG0745-COG2203-COG2203-COG5001','COG2202-COG3829-COG5001','COG0715-COG3706','COG2199-COG2200-COG3829-COG3829','COG2199-COG2200-COG4753','COG3300-COG5001','COG3275-COG5001','COG3437-COG3829-COG4191-COG5001','COG2202-COG2203-COG3706','COG2199-COG2200-COG3829-COG3852-COG5002','COG0784-COG3290-COG3706','COG2202-COG5001-COG5002','COG2199-COG2200-COG3710','COG0683-COG3290-COG5001','COG2703-COG3706','COG2703-COG5001','COG4252-COG5001','COG0784-COG2199-COG2200','COG3614-COG3706','COG3706-COG5082','COG3706-COG3850','COG2203-COG2461-COG3706','COG2199-COG2202','COG2203-COG2203-COG3706','COG3447-COG4191-COG5001','COG2202-COG3221-COG5001-COG5002','COG2199-COG2703-COG4191','COG2203-COG3706-COG5000','COG0517-COG0517-COG3829-COG5001','COG2199-COG2202-COG2203','COG1086-COG2199-COG5000','COG3706-COG3850-COG4564','COG2204-COG3829-COG5001')

#EAL
COGselect_grep = "EAL"

COGselect = c('COG3437-COG5001','COG3437-COG3829-COG5001','COG5001','COG0745-COG3290-COG5001','COG4191-COG5001','COG4564-COG5001','COG2202-COG2202-COG5001','COG3706-COG5001','COG2203-COG5001','COG3434','COG0745-COG5001','COG0784-COG0784-COG2198-COG2205-COG3437-COG5001','COG2202-COG5001','COG4191-COG5000-COG5001','COG2905-COG2905-COG5001','COG2202-COG3614-COG5001','COG5000-COG5001-COG5002','COG3852-COG5001','COG2202-COG4251-COG5001','COG2202-COG5000-COG5001','COG3447-COG3829-COG5001','COG2202-COG3706-COG5001','COG5001-COG5002','COG3850-COG5001','COG2202-COG3284-COG5001','COG3829-COG5001','COG3829-COG5001-COG5002','COG1196-COG1352-COG2201-COG3829-COG3829-COG5001','COG0745-COG2202-COG5001-COG5002','COG1352-COG2201-COG3829-COG4717-COG5001','COG2202-COG2202-COG2203-COG5001','COG2203-COG5000-COG5001','COG3452-COG5001','COG2203-COG2203-COG5001','COG2203-COG2203-COG2203-COG2203-COG2203-COG5001','COG1716-COG5001','COG2205-COG3829-COG5001','COG5000-COG5001','COG2202-COG2203-COG3290-COG5001','COG2202-COG2203-COG5001','COG0457-COG2202-COG5001','COG2202-COG2461-COG5001','COG0004-COG2202-COG3829-COG5000-COG5001','COG3284-COG5001','COG0784-COG2203-COG2203-COG5001','COG2203-COG4191-COG5001','COG2203-COG3290-COG5001','COG2202-COG3279-COG5001','COG0784-COG2198-COG2202-COG5001','COG2203-COG3829-COG5001','COG2202-COG2202-COG2203-COG2203-COG2203-COG5001','COG4943','COG3614-COG5001','COG0745-COG2203-COG2203-COG5001','COG2202-COG3829-COG5001','COG3300-COG5001','COG3275-COG5001','COG3437-COG3829-COG4191-COG5001','COG2202-COG5001-COG5002','COG0683-COG3290-COG5001','COG2703-COG5001','COG4252-COG5001','COG3447-COG4191-COG5001','COG2202-COG3221-COG5001-COG5002','COG0517-COG0517-COG3829-COG5001','COG2204-COG3829-COG5001')
#GGDEF, EAL, HD
COGselect_grep = "GGDEF, EAL, HD"

COGselect = c('COG3437-COG5001','COG3437-COG3829-COG5001','COG3437','COG0784-COG0784-COG2198-COG2205-COG3437-COG5001','COG3437-COG3437','COG0745-COG2114-COG2203-COG2203-COG3437-COG5002','COG0784-COG3437-COG4191','COG2203-COG2205-COG3437','COG2203-COG3437','COG2206-COG3437','COG2114-COG2203-COG2203-COG3437','COG0834-COG2205-COG3437','COG2206','COG3437-COG3829-COG4191-COG5001','COG3437-COG4191','COG2804-COG3437','COG2199-COG2200','COG3437-COG5001','COG3437-COG3829-COG5001','COG2200-COG3706','COG0517-COG2199-COG2200','COG2199-COG2203','COG3437','COG3434','COG2199','COG2199-COG2202-COG2203-COG3852-COG5002','COG0784-COG0784-COG2198-COG2205-COG3437-COG5001','COG2199-COG2200-COG3287','COG2199-COG2203-COG3829','COG2200','COG2199-COG2703','COG2199-COG2200-COG2204','COG2199-COG4191','COG2199-COG2905-COG2905','COG3437-COG3437','COG0745-COG2114-COG2203-COG2203-COG3437-COG5002','COG2199-COG2200-COG2203','COG2200-COG2203','COG0784-COG3437-COG4191','COG2199-COG2200-COG2203-COG3290','COG2199-COG2200-COG3279','COG2203-COG2205-COG3437','COG2203-COG3437','COG2199-COG3300','COG2199-COG2200-COG3290','COG0784-COG2200-COG4566','COG2200-COG2204-COG3706','COG1835-COG2199','COG2199-COG3829','COG2206-COG3437','COG2114-COG2203-COG2203-COG3437','COG2199-COG2200-COG3850','COG2199-COG3447','COG0697-COG0784-COG2198-COG2199-COG2200-COG2205','COG0834-COG2205-COG3437','COG2199-COG3221-COG5002','COG2206','COG4943','COG2199-COG2200-COG3829-COG3829','COG2199-COG2200-COG4753','COG0745-COG2200','COG3437-COG3829-COG4191-COG5001','COG3437-COG4191','COG2804-COG3437','COG2199-COG2200-COG3829-COG3852-COG5002','COG2199-COG2200-COG3710','COG0784-COG2199-COG2200','COG2199-COG2202','COG2199-COG2703-COG4191','COG2199-COG2202-COG2203','COG1086-COG2199-COG5000','COG3437-COG5001','COG3437-COG3829-COG5001','COG5001','COG0745-COG3290-COG5001','COG4191-COG5001','COG4564-COG5001','COG2202-COG2202-COG5001','COG3706-COG5001','COG2203-COG5001','COG3434','COG0745-COG5001','COG0784-COG0784-COG2198-COG2205-COG3437-COG5001','COG2202-COG5001','COG4191-COG5000-COG5001','COG2905-COG2905-COG5001','COG2202-COG3614-COG5001','COG5000-COG5001-COG5002','COG3852-COG5001','COG2202-COG4251-COG5001','COG2202-COG5000-COG5001','COG3447-COG3829-COG5001','COG2202-COG3706-COG5001','COG5001-COG5002','COG3850-COG5001','COG2202-COG3284-COG5001','COG3829-COG5001','COG3829-COG5001-COG5002','COG1196-COG1352-COG2201-COG3829-COG3829-COG5001','COG0745-COG2202-COG5001-COG5002','COG1352-COG2201-COG3829-COG4717-COG5001','COG2202-COG2202-COG2203-COG5001','COG2203-COG5000-COG5001','COG3452-COG5001','COG2203-COG2203-COG5001','COG2203-COG2203-COG2203-COG2203-COG2203-COG5001','COG1716-COG5001','COG2205-COG3829-COG5001','COG5000-COG5001','COG2202-COG2203-COG3290-COG5001','COG2202-COG2203-COG5001','COG0457-COG2202-COG5001','COG2202-COG2461-COG5001','COG0004-COG2202-COG3829-COG5000-COG5001','COG3284-COG5001','COG0784-COG2203-COG2203-COG5001','COG2203-COG4191-COG5001','COG2203-COG3290-COG5001','COG2202-COG3279-COG5001','COG0784-COG2198-COG2202-COG5001','COG2203-COG3829-COG5001','COG2202-COG2202-COG2203-COG2203-COG2203-COG5001','COG4943','COG3614-COG5001','COG0745-COG2203-COG2203-COG5001','COG2202-COG3829-COG5001','COG3300-COG5001','COG3275-COG5001','COG3437-COG3829-COG4191-COG5001','COG2202-COG5001-COG5002','COG0683-COG3290-COG5001','COG2703-COG5001','COG4252-COG5001','COG3447-COG4191-COG5001','COG2202-COG3221-COG5001-COG5002','COG0517-COG0517-COG3829-COG5001','COG2204-COG3829-COG5001','COG0784-COG2198-COG2202-COG2205-COG3290-COG3706-COG3706-COG3829','COG2199-COG2200','COG3437-COG5001','COG0745-COG2198-COG3706','COG3437-COG3829-COG5001','COG3706-COG5002','COG5001','COG3706','COG0745-COG3290-COG5001','COG2200-COG3706','COG2202-COG2202-COG3706','COG4191-COG5001','COG4564-COG5001','COG0517-COG2199-COG2200','COG2199-COG2203','COG2202-COG2202-COG5001','COG3706-COG5001','COG2203-COG5001','COG2199','COG2199-COG2202-COG2203-COG3852-COG5002','COG0745-COG5001','COG0784-COG0784-COG2198-COG2205-COG3437-COG5001','COG2202-COG5001','COG4191-COG5000-COG5001','COG2905-COG2905-COG5001','COG2199-COG2200-COG3287','COG2199-COG2203-COG3829','COG3706-COG3829','COG2202-COG3614-COG5001','COG5000-COG5001-COG5002','COG2199-COG2703','COG0697-COG3706','COG2199-COG2200-COG2204','COG3322-COG3706','COG3852-COG5001','COG2199-COG4191','COG3706-COG5000','COG2202-COG4251-COG5001','COG2202-COG3706','COG2202-COG5000-COG5001','COG3447-COG3829-COG5001','COG2203-COG3300-COG3706','COG2199-COG2905-COG2905','COG2202-COG3706-COG5001','COG5001-COG5002','COG3850-COG5001','COG2202-COG3284-COG5001','COG2199-COG2200-COG2203','COG3829-COG5001','COG3829-COG5001-COG5002','COG1196-COG1352-COG2201-COG3829-COG3829-COG5001','COG0745-COG2202-COG5001-COG5002','COG1352-COG2201-COG3829-COG4717-COG5001','COG2202-COG2202-COG2203-COG5001','COG2203-COG5000-COG5001','COG2203-COG3706','COG3452-COG5001','COG2203-COG2203-COG5001','COG3706-COG4251','COG2203-COG2203-COG2203-COG2203-COG2203-COG5001','COG1716-COG5001','COG2205-COG3829-COG5001','COG5000-COG5001','COG2199-COG2200-COG2203-COG3290','COG0834-COG0834-COG0834-COG3706','COG2202-COG2203-COG3290-COG5001','COG2202-COG2203-COG5001','COG0834-COG0834-COG3706','COG2199-COG2200-COG3279','COG0457-COG2202-COG5001','COG2202-COG2461-COG5001','COG0004-COG2202-COG3829-COG5000-COG5001','COG3284-COG5001','COG0784-COG2203-COG2203-COG5001','COG2199-COG3300','COG2203-COG4191-COG5001','COG2199-COG2200-COG3290','COG3706-COG3852','COG2200-COG2204-COG3706','COG2203-COG3290-COG5001','COG1835-COG2199','COG2199-COG3829','COG0457-COG0515-COG2203-COG2203-COG3706-COG3899','COG0835-COG3706','COG2198-COG3706','COG2202-COG3279-COG5001','COG0784-COG2198-COG2202-COG5001','COG0745-COG3706','COG2199-COG2200-COG3850','COG2199-COG3447','COG3264-COG3706','COG0697-COG0784-COG2198-COG2199-COG2200-COG2205','COG3706-COG4191','COG0715-COG2703-COG3706','COG0784-COG2202-COG2202-COG3706','COG2208-COG3706','COG3447-COG3706','COG2199-COG3221-COG5002','COG2203-COG3829-COG5001','COG2202-COG2202-COG2203-COG2203-COG2203-COG5001','COG0834-COG3706','COG3614-COG5001','COG3706-COG4963','COG0745-COG2203-COG2203-COG5001','COG2202-COG3829-COG5001','COG0715-COG3706','COG2199-COG2200-COG3829-COG3829','COG2199-COG2200-COG4753','COG3300-COG5001','COG3275-COG5001','COG3437-COG3829-COG4191-COG5001','COG2202-COG2203-COG3706','COG2199-COG2200-COG3829-COG3852-COG5002','COG0784-COG3290-COG3706','COG2202-COG5001-COG5002','COG2199-COG2200-COG3710','COG0683-COG3290-COG5001','COG2703-COG3706','COG2703-COG5001','COG4252-COG5001','COG0784-COG2199-COG2200','COG3614-COG3706','COG3706-COG5082','COG3706-COG3850','COG2203-COG2461-COG3706','COG2199-COG2202','COG2203-COG2203-COG3706','COG3447-COG4191-COG5001','COG2202-COG3221-COG5001-COG5002','COG2199-COG2703-COG4191','COG2203-COG3706-COG5000','COG0517-COG0517-COG3829-COG5001','COG2199-COG2202-COG2203','COG1086-COG2199-COG5000','COG3706-COG3850-COG4564','COG2204-COG3829-COG5001')


#arginine
COGselect_grep = "arginine"

COGselect = c('COG0810','COG0810-COG3291','COG0810-COG4219')

#ferredoxin and complex1
COGselect_grep = "ferredoxin and complex1"

COGselect = c('COG0543-COG0633','COG0348','COG1018','COG0348-COG3901','COG4656-COG4656','COG2146','COG0633','COG1014-COG4231','COG2878','COG1145','COG0543-COG1018','COG4659','COG0243-COG2906','COG0543-COG1018-COG3239','COG2906','COG1018-COG3576-COG3576','COG1018-COG1018','COG1018-COG3000','COG1145-COG1149','COG1145-COG1148','COG1139','COG1017-COG1018','COG0247-COG0277-COG1145','COG4660','COG4658','COG4657','COG4656','COG0493-COG1144','COG0674-COG1013-COG1014-COG1143','COG0674-COG1013-COG1014-COG1146','COG0644-COG2440','COG1018-COG2124','COG1842-COG4656','COG1196-COG4656','COG1141','COG1018-COG3576','COG0838','COG0377','COG0852','COG0649','COG1905','COG1894','COG1034','COG1005','COG1143','COG0839','COG0713','COG1009','COG1008','COG1007')

#Complex1
COGselect_grep = "complex1"

COGselect = c('COG0838','COG0377','COG0852','COG0649','COG1905','COG1894','COG1034','COG1005','COG1143','COG0839','COG0713','COG1009','COG1008','COG1007')

#beta-oxidation
COGselect_grep = "beta-oxidation"

COGselect = c('COG1960','COG1024','COG1250','COG0183')


#grep NAD
COGselect_grep = "NAD"

COGselect = c('COG1034','COG1012','COG0702','COG1028','COG0543-COG0633','COG0713','COG1748','COG1231-COG1902','COG1902','COG0604','COG1007','COG0846','COG2902','COG1018','COG3935','COG0272','COG1251','COG4221','COG4656-COG4656','COG0431','COG0446','COG2878','COG0377','COG1009','COG0543','COG2130','COG0543-COG1018','COG1008','COG0780-COG2904','COG4659','COG0243-COG1251','COG0569','COG1090','COG1252','COG0272-COG0507','COG1182','COG1726','COG2907','COG0838','COG2070','COG0569-COG0569','COG0543-COG1018-COG3239','COG0839','COG2249','COG0404-COG0446','COG3367','COG1018-COG3576-COG3576','COG0236-COG0300-COG0604-COG3321','COG1018-COG1018','COG3380','COG3236','COG0446-COG1902','COG1018-COG3000','COG2871','COG3288','COG1005','COG0272-COG0847','COG0471-COG0471-COG0569','COG0649-COG0852','COG0236-COG0236-COG0236-COG0236-COG0236-COG0236-COG0300-COG0500-COG1024-COG1028-COG3321-COG3321-COG3321-COG5644','COG1017-COG1018','COG1853','COG0493','COG1146','COG0471-COG0471-COG0569-COG3273','COG1009-COG2111','COG0171','COG1805-COG1805','COG0623','COG4660','COG4658','COG2869','COG1347','COG4657','COG2375','COG1073-COG4221','COG1012-COG2030','COG1859','COG0446-COG1773','COG1251-COG1773','COG1028-COG3320','COG4656','COG1143','COG0061','COG0493-COG1144','COG1846-COG1853','COG1282','COG0702-COG0702','COG0062-COG0063','COG2816','COG2209','COG0846-COG1611','COG0446-COG3453','COG0025-COG0569','COG1905','COG2344','COG0674-COG1013-COG1014-COG1143','COG0674-COG1013-COG1014-COG1146','COG1009-COG2111-COG2111','COG1894-COG1905','COG1018-COG2124','COG0479-COG0543','COG0709-COG1252','COG0382-COG2249','COG1894','COG1842-COG4656','COG1196-COG4656','COG0247-COG1146','COG1018-COG3576','COG0171-COG0388','COG0543-COG0633','COG0348','COG1018','COG0348-COG3901','COG4656-COG4656','COG2146','COG0633','COG1014-COG4231','COG2878','COG1145','COG0543-COG1018','COG4659','COG0243-COG2906','COG0543-COG1018-COG3239','COG2906','COG1018-COG3576-COG3576','COG1018-COG1018','COG1018-COG3000','COG1145-COG1149','COG1145-COG1148','COG1139','COG1017-COG1018','COG0247-COG0277-COG1145','COG4660','COG4658','COG4657','COG4656','COG0493-COG1144','COG0674-COG1013-COG1014-COG1143','COG0674-COG1013-COG1014-COG1146','COG0644-COG2440','COG1018-COG2124','COG1842-COG4656','COG1196-COG4656','COG1141','COG1018-COG3576','COG0838','COG0377','COG0852','COG0649','COG1905','COG1894','COG1034','COG1005','COG1143','COG0839','COG0713','COG1009','COG1008','COG1007','COG1960','COG1024','COG1250','COG0183')

#regulatory
COGselect = c('COG1309','COG0640','COG0583','COG3182','COG1737','COG1585','COG2197','COG0745','COG0745-COG1639','COG0784-COG2198-COG2202-COG2205-COG3290-COG3706-COG3706-COG3829','COG2865','COG1802','COG0745-COG2197','COG4191','COG1401','COG3026','COG0789','COG4977','COG3835','COG2201','COG3437-COG5001','COG0745-COG2198-COG3706','COG1476','COG3437-COG3829-COG5001','COG3706-COG5002','COG1522','COG3829-COG4251','COG5631','COG3705','COG1678','COG2204','COG4467','COG2771','COG1396-COG2856','COG3706','COG1396','COG2001','COG2188','COG1349','COG0745-COG3290-COG5001','COG2200-COG3706','COG2202-COG2202-COG3706','COG4191-COG5001','COG3311','COG0840-COG3829','COG0840-COG3829-COG5000','COG2172','COG2203-COG3300-COG4251','COG0628','COG2378','COG0741','COG1733','COG1414','COG0745-COG5002','COG5000','COG5340','COG0664','COG4191-COG4191','COG3279','COG3706-COG5001','COG3905','COG0348-COG3901','COG1401-COG4127','COG3437','COG2202-COG2203-COG2203-COG2205-COG3829','COG3829','COG0745-COG5001','COG0784-COG0784-COG2198-COG2205-COG3437-COG5001','COG4191-COG5000-COG5001','COG1327','COG3604','COG1609','COG0330','COG2202-COG4191','COG0745-COG0745-COG2202-COG2205','COG2199-COG2203-COG3829','COG2172-COG2208-COG2972','COG5662','COG0347','COG3706-COG3829','COG0643-COG0745-COG2198-COG2198-COG2198-COG2198-COG2198-COG2198','COG0661','COG3850-COG4251','COG5000-COG5001-COG5002','COG1739','COG2390','COG2186','COG1551','COG0745-COG0784-COG2198-COG5002','COG1937','COG0596-COG2197','COG3076','COG0457-COG1729-COG1729-COG1729-COG2956','COG1366','COG0745-COG2202-COG2205-COG4191','COG1222-COG1222','COG1846','COG3609','COG0745-COG1457-COG2205','COG0697-COG3706','COG2110','COG2909','COG1396-COG1917','COG0323-COG4191','COG2199-COG2200-COG2204','COG3290','COG3322-COG3706','COG3118','COG2199-COG4191','COG0664-COG0745','COG3655','COG0834-COG4191','COG3706-COG5000','COG2202-COG4251-COG5001','COG0784-COG2205-COG4251','COG4282','COG0840-COG3290','COG0784-COG4251','COG2202-COG3706','COG0784-COG2202-COG2202-COG2203-COG4191','COG2202-COG5000-COG5001','COG3447-COG3829-COG5001','COG0270-COG2944','COG2203-COG3300-COG3706','COG2203-COG2208','COG1917-COG4977','COG2202-COG3706-COG5001','COG3437-COG3437','COG4268','COG1281','COG0745-COG2114-COG2203-COG2203-COG3437-COG5002','COG3057','COG2202-COG3284-COG5001','COG2005','COG0378','COG1522-COG1522','COG0789-COG2226','COG0789-COG2230','COG3829-COG5001','COG3707','COG3829-COG5001-COG5002','COG4197','COG1695','COG4933','COG1479-COG2865-COG3472','COG0369-COG3182','COG3614-COG4251','COG2771-COG3292','COG2915','COG3712','COG1196-COG1352-COG2201-COG3829-COG3829-COG5001','COG2203-COG4191','COG0745-COG2202-COG5001-COG5002','COG1352-COG2201-COG3829-COG4717-COG5001','COG4566','COG1167','COG1510','COG3016','COG3070','COG2203-COG5000-COG5001','COG3086','COG2203-COG3706','COG3180','COG2002','COG0697-COG0745-COG4191-COG5002','COG3226','COG3423','COG2197-COG2267','COG3228','COG0664-COG1639','COG2808','COG3706-COG4251','COG2944','COG2205-COG3829-COG5001','COG3284','COG1959','COG0784-COG3437-COG4191','COG5000-COG5001','COG2197-COG2202-COG2202-COG4191-COG4936','COG1396-COG2932','COG2199-COG2200-COG2203-COG3290','COG2018','COG4753','COG0834-COG0834-COG0834-COG3706','COG2202-COG2203-COG3290-COG5001','COG0741-COG0791','COG0704','COG5450','COG5566','COG0735','COG0834-COG0834-COG3706','COG1396-COG3837','COG1719-COG3829','COG5602','COG3327','COG0784-COG2205-COG5000','COG2199-COG2200-COG3279','COG0642-COG0745-COG2202-COG2202-COG2203-COG3829-COG5000','COG0784-COG2202-COG2203-COG5000-COG5002','COG0004-COG2202-COG3829-COG5000-COG5001','COG0318-COG1309','COG3284-COG5001','COG3620','COG2203-COG2205-COG3437','COG2203-COG3437','COG1846-COG3177','COG0741-COG3941','COG2203-COG4191-COG5001','COG2199-COG2200-COG3290','COG2197-COG2865','COG3706-COG3852','COG1476-COG2856','COG4565-COG4584','COG2771-COG4584','COG1806','COG0784-COG2200-COG4566','COG2200-COG2204-COG3706','COG3357','COG0745-COG2208','COG0745-COG2909','COG2508','COG0745-COG2207','COG1776-COG2201','COG0684','COG3703','COG3487','COG0784-COG2202-COG2202-COG4191','COG2203-COG3290-COG5001','COG1598-COG3609','COG0745-COG3447-COG3614-COG5002','COG3682','COG0217','COG1716-COG2204','COG0591-COG4191','COG2202-COG2202-COG2202-COG3829-COG4251','COG3937','COG2199-COG3829','COG0457-COG0515-COG2203-COG2203-COG3706-COG3899','COG3722','COG0471-COG0471-COG0490','COG0835-COG3706','COG5007','COG4935','COG2198-COG3706','COG1366-COG3852','COG2206-COG3437','COG4567','COG0784-COG4191','COG2202-COG3279-COG5001','COG0864-COG1598','COG2207-COG3708','COG2114-COG2203-COG2203-COG3437','COG0745-COG2197-COG2202-COG2202-COG3829-COG4251','COG3646','COG0745-COG3706','COG2204-COG2207','COG3264-COG3706','COG0745-COG2205','COG0784-COG0784-COG2198-COG3850-COG4251','COG0834-COG2205-COG3437','COG3022','COG2964','COG3706-COG4191','COG0745-COG0784-COG2198-COG2202-COG2202-COG3852-COG4251-COG5000','COG0784-COG0784-COG2198-COG2202-COG2202-COG2205-COG3829-COG5002','COG0715-COG2703-COG3706','COG2976','COG0784-COG2202-COG2203-COG2205-COG3614-COG4564-COG5000','COG0784-COG2202-COG2202-COG3706','COG0394-COG0640','COG2208-COG3706','COG4191-COG4192','COG3160','COG3447-COG3706','COG1846-COG1853','COG0784-COG2197-COG2205','COG0810-COG4219','COG2852-COG5340','COG2203-COG3829-COG5001','COG2202-COG2202-COG2202-COG2202-COG2203-COG2203-COG3290-COG3290-COG3290-COG4251','COG2137','COG2203-COG3279','COG0789-COG1192','COG2319-COG2319-COG2319-COG2909','COG0745-COG4977','COG0149-COG3707','COG1664-COG5662','COG0750','COG5394','COG1910-COG2005','COG2522','COG0640-COG4976','COG0515-COG2203-COG3899-COG4191','COG0834-COG3706','COG4251','COG1196-COG1401','COG0607-COG0640','COG2956','COG3250-COG4189','COG4572','COG4650','COG3706-COG4963','COG0741-COG1388-COG1388-COG1388','COG0745-COG2203-COG2203-COG5001','COG2202-COG3829-COG5001','COG0745-COG0784-COG2198-COG4251','COG0715-COG3706','COG0640-COG1708','COG2344','COG3073','COG4565','COG2199-COG2200-COG3829-COG3829','COG2199-COG2200-COG4753','COG0784-COG2197-COG2198-COG2205','COG2203-COG3829','COG2747','COG0745-COG2200','COG0864','COG2202-COG4251','COG2202-COG3829-COG4251-COG5002','COG3322-COG4192-COG4251','COG0475-COG0490-COG1226','COG3437-COG3829-COG4191-COG5001','COG1395','COG0457-COG2909','COG0789-COG2442','COG2921','COG0666-COG0666-COG0666-COG4282','COG3437-COG4191','COG2804-COG3437','COG2202-COG2203-COG3706','COG2199-COG2200-COG3829-COG3852-COG5002','COG0784-COG3290-COG3706','COG2208','COG0683-COG3290-COG5001','COG2703-COG3706','COG3143','COG3614-COG3706','COG3706-COG5082','COG3829-COG4585','COG3706-COG3850','COG1401-COG2256','COG0470-COG1401','COG1923','COG2203-COG2461-COG3706','COG0642-COG5000','COG1802-COG1802','COG2203-COG2203-COG3706','COG1983','COG3447-COG4191-COG5001','COG1396-COG1680-COG3296','COG0745-COG2205-COG5278','COG2186-COG3618','COG2199-COG2703-COG4191','COG2203-COG3706-COG5000','COG0517-COG0517-COG3829-COG5001','COG0799','COG0758-COG1846','COG0643-COG0643-COG0745-COG2198','COG1086-COG2199-COG5000','COG1196-COG1352-COG2201','COG0784-COG2202-COG4191','COG2202-COG2202-COG4191-COG4251','COG3706-COG3850-COG4564','COG0457-COG2956','COG1352-COG2201-COG3352','COG1352-COG2201-COG2461','COG2204-COG3829-COG5001','COG3614-COG3829-COG5002')

#resp-ferm
COGselect_grep = "resp-ferm"
COGselect = c('COG1073','COG1073-COG2267','COG1073-COG4221','COG0508-COG1073','COG1073-COG1765')

#nitrous
COGselect_grep = "denitrification (key)"
COGselect = c('COG4263','COG5013',"COG2146",'COG3256')

#nitrous
COGselect_grep = "NADH-ubiquinone"
COGselect = c('COG0713','COG1007','COG1009','COG1008','COG0838','COG0839','COG1005')
#nitrous
COGselect_grep = "formate"
COGselect = c('COG0651','COG0027','COG1526','COG2864','COG2116','COG1143','COG1180','COG0457-COG1413-COG3303','COG0674-COG1013-COG1014-COG1143','COG3058')
COGselect_grep = "sodium_symport"
COGselect = c('COG2111','COG0651','COG2169','COG1115','COG3135','COG0025-COG1226','COG1283','COG4656-COG4656','COG0591-COG0591','COG0591','COG0350-COG2169','COG3630','COG2878','COG1009','COG0591-COG2205','COG4659','COG0534','COG0025','COG1320','COG2610','COG1726','COG0733','COG0385','COG1863','COG2212','COG4145','COG3263','COG1301','COG2871','COG1006','COG4146','COG1055','COG2851','COG0350','COG1009-COG2111','COG1805-COG1805','COG4660','COG0591-COG4191','COG4658','COG2869','COG1347','COG4147','COG4657','COG0122-COG2169','COG1883','COG4656','COG0786','COG2935','COG2209','COG3633','COG0025-COG0569','COG2211','COG1009-COG2111-COG2111','COG1757','COG0591-COG0784-COG2205','COG3067','COG1842-COG4656','COG1196-COG4656','COG0025-COG2905','COG3004')
COGselect_grep = "cold-heat"
COGselect = c('COG3187','COG1278','COG1413','COG3187-COG3650-COG3895','COG3187-COG3187','COG0457-COG1413','COG3015-COG3187','COG1278-COG3326','COG1413-COG2010-COG2133','COG3126-COG3187','COG0576','COG0457-COG1413-COG3303','COG1409-COG1413-COG1413-COG5635','COG0323-COG1413','COG0559-COG1413','COG1188','COG1413-COG1413-COG4249-COG5635')
COGselect_grep = "lactate"
COGselect = c('COG0028','COG1052','COG0039','COG1620','COG1139','COG2055','COG0440','COG1304','COG3978','COG0391','COG1920','COG1556')
COGselect_grep = "glucose"
COGselect = c('COG2133','COG0357','COG0166','COG2133-COG3794','COG2010-COG4993','COG2246','COG4993','COG1087','COG2133-COG3291-COG3828-COG4654','COG2133-COG3055','COG0364','COG1413-COG2010-COG2133','COG1210','COG1004','COG1209','COG0448','COG1215-COG2246','COG1088')
COGselect_grep = "NADPH"
COGselect = c('COG0604','COG0446','COG2130','COG0780-COG2904','COG2249','COG0404-COG0446','COG0236-COG0300-COG0604-COG3321','COG0446-COG1902','COG0493','COG2375','COG0446-COG1773','COG0493-COG1144','COG0446-COG3453','COG0382-COG2249')
COGselect_grep = "FAdesaturase"
COGselect = c('COG3239','COG3000','COG0543-COG1018-COG3239','COG1018-COG3000','COG1398','COG1398-COG3464')
COGselect_grep = "glutam"
COGselect = c('COG4948','COG0367','COG0548','COG2201','COG0334','COG0518','COG3828','COG2902','COG0796','COG1305','COG0069','COG3931','COG3653','COG2071','COG1794','COG0405','COG2843','COG2988','COG0285','COG4865','COG2170','COG0008','COG1871','COG2133-COG3291-COG3828-COG4654','COG0174','COG3672','COG1196-COG1352-COG2201-COG3829-COG3829-COG5001','COG0118','COG2105','COG1352-COG2201-COG3829-COG4717-COG5001','COG1478','COG0046-COG0047','COG2918','COG1305-COG4196','COG0034','COG0518-COG0519','COG1680-COG3653','COG0067-COG0069-COG0070','COG1800','COG0002','COG0697-COG1305','COG0014','COG0493','COG1776-COG2201','COG1305-COG1372-COG3209-COG3209-COG3209','COG0001-COG1020-COG3321','COG0001-COG0236-COG1020-COG3321','COG0076','COG0121','COG0373','COG0069-COG3369-COG3369','COG0493-COG1144','COG0001','COG0786','COG1391','COG0747-COG1305','COG0548-COG1246','COG0771','COG0263','COG1364','COG2066','COG0034-COG0367','COG1196-COG1352-COG2201','COG1352-COG2201-COG3352','COG1352-COG2201-COG2461')

#nitrous
COGselect_grep = "oxidoreductase"
COGselect = c('COG1034','COG0654','COG0713','COG2141','COG0247-COG0277','COG0604','COG1007','COG0667','COG4656-COG4656','COG1014-COG4231','COG0446','COG2878','COG3011','COG0377','COG1009','COG1008','COG4659','COG3383','COG0635','COG1726','COG0838','COG3573','COG2070','COG0839','COG0404-COG0446','COG0236-COG0300-COG0604-COG3321','COG0243-COG3383','COG3380','COG0446-COG1902','COG2871','COG0247','COG1005','COG0649-COG0852','COG1853','COG4312','COG0493','COG1009-COG2111','COG1805-COG1805','COG0247-COG0277-COG1145','COG0247-COG0277-COG0479','COG4660','COG4658','COG2869','COG1347','COG0667-COG1861','COG3560','COG4657','COG0446-COG1773','COG4656','COG1143','COG0493-COG1144','COG1846-COG1853','COG2209','COG0446-COG3453','COG1905','COG0674-COG1013-COG1014-COG1143','COG0674-COG1013-COG1014-COG1146','COG1453','COG4989','COG1009-COG2111-COG2111','COG1894-COG1905','COG1894','COG1842-COG4656','COG1196-COG4656','COG0247-COG1146','COG0247-COG0477')

#nitrous
COGselect_grep = "FAD"
COGselect = c('COG0654','COG0247-COG0277','COG0175','COG1319','COG0824','COG0175-COG3415','COG2186','COG2509','COG0277','COG1252','COG2907','COG3380','COG0196','COG0247-COG0277-COG1145','COG0247-COG0277-COG0479','COG2375','COG4630','COG0709-COG1252','COG0196-COG3172')
COGselect_grep = "flagellum"
COGselect = c('COG3144','COG1360','COG2063','COG1345-COG1345','COG1338','COG1886','COG1419','COG0455','COG1256','COG1298','COG1536','COG4787','COG1815','COG2882','COG1157','COG1291','COG1344-COG1344','COG1344','COG1377','COG1766','COG1684','COG1580','COG3190','COG1843-COG3656','COG1516','COG1677','COG1558','COG1868','COG1344-COG1409','COG1317','COG2747','COG1749','COG3418','COG1706','COG1705-COG3951','COG4786','COG1987','COG1261')
COGselect_grep = "T6SS"
COGselect = c('COG3520','COG3157','COG3501','COG3517','COG3519','COG3522','COG3523','COG3521','COG3516','COG3515-COG5283-COG5412','COG3455','COG3515','COG3518','COG1716-COG3456')
COGselect_grep = "TonB"
COGselect = c('COG0810','COG0810-COG3291','COG0810-COG4219')
COGselect_grep = "hydrolase"
COGselect = c('COG0491','COG2333','COG3179-COG3409','COG0388','COG2267','COG1574','COG3409','COG0791','COG4942','COG4757','COG2159','COG0402','COG0412','COG3571','COG1835','COG1451','COG5533','COG3409-COG4322','COG0739','COG1524-COG3568','COG2049','COG3931','COG2071','COG3687','COG0658-COG2333','COG0499','COG1075','COG3453','COG1228','COG0427-COG0427','COG3621','COG1397','COG1607','COG0500-COG0537','COG1073-COG2267','COG5621','COG0491-COG0607','COG1078','COG3194','COG0494','COG1402','COG0788','COG1408','COG3964','COG4927','COG0327','COG2351','COG1984','COG0807','COG2197-COG2267','COG0429','COG0327-COG0327-COG3323','COG0427','COG0232','COG0302','COG0193','COG0741-COG0791','COG0388-COG1247','COG0439-COG0511-COG1984-COG2049','COG1397-COG1397-COG1397','COG5520','COG0739-COG3061','COG2015','COG1380','COG0491-COG0607-COG0607','COG1896','COG3409-COG3772','COG1346','COG0537','COG0627','COG1469','COG3356','COG2503-COG3568','COG0190','COG4339','COG1835-COG2199','COG2248','COG0044','COG1388-COG4942','COG4341','COG0388-COG3344','COG3115-COG3409','COG0797-COG3087','COG0139','COG1387-COG1796','COG0352-COG0494','COG3568','COG4188','COG2951-COG3409','COG0388-COG3153','COG0739-COG1388','COG1988','COG2429','COG0248','COG2175-COG4341','COG0446-COG3453','COG3970','COG1835-COG2755','COG0599-COG2267','COG0599-COG0599-COG1075-COG1917','COG4225','COG1554','COG1621','COG0138','COG3545','COG3618','COG2819','COG0797','COG1835-COG4763','COG3724','COG2152','COG1480','COG1501','COG0140','COG0457-COG2819','COG1705-COG3951','COG0108-COG0807','COG2945','COG2186-COG3618','COG0210-COG1112-COG1502-COG4942','COG1957','COG0171-COG0388','COG0491-COG0494')
########
COGselect_grep = "denitrification (key)"
COGselect = c('COG4263','COG5013',"COG2146",'COG3256')
COGselect = head(finRes2$COG_FUNCTION_ACC,50)
DNF = ifelse(rowSums(mat)>2,"denitrifyer","non-denitrifyer")
DNF = data.frame(DNF)
DNF$genome = rownames(DNF)
dummy = data.frame(genome = levels(ANVIO$genome))
DNF2 = right_join(DNF,dummy, by=c('genome'))
DNF2$DNF[is.na(DNF2$DNF)]="non-denitrifyer"
colnames(DNF2) = c('DNF',"Bin_Id") 

###########





#mat = t(table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% COGselect)[,c('COG_FUNCTION_ACC','genome_name')])))

#mat = t(table(droplevels(ANVIO[grepl('regula', ANVIO$COG_FUNCTION),c('COG_FUNCTION','genome_name')])))

mat = t(table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% COGselect)[,c('COG_FUNCTION','genome_name')])))

###########
#CHECK AND RESOLVE ISSUE
###########
mat = t(table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% COGselect)[,c('COG_FUNCTION','genome_name')])))



dmat = data.frame(mat)
dmat = reshape(dmat, idvar = "genome_name", timevar = "COG_FUNCTION", direction = "wide")

colnames(dmat) = gsub("Freq.",'',colnames(dmat))
rownames(dmat) = dmat $genome_name
dmat = dmat[,2:length(colnames(dmat))]
dmat$genome = rownames(dmat)
matFull = right_join(dmat,dummy, by=c('genome'))
matFull[is.na(matFull)]=0
rownames(matFull) = matFull$genome
matFull = matFull[,1:length(colnames(matFull))-1]
mat = matFull
###########
###########
###########

hhist1 = ggplot(data.frame(rowSums(mat)),aes(x=rowSums.mat.)) + geom_histogram( colour="black", fill="grey") + geom_vline(aes(xintercept=mean(rowSums.mat.)),color="grey", linetype="dashed", size=0.75) + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + xlab(paste('nr of ', COGselect_grep,sep='')) + theme_classic()

transposases = data.frame(rowSums(mat))
transposases = data.frame(rowSums(mat), Bin_Id = rownames(transposases))
Tps = right_join(transposases, CheckM2, by=c('Bin_Id'))

ANVIO_cat2=ANVIO_cat
names(ANVIO_cat2) = c("Bin_Id","group","lifestyle","SS","sourceClass","lat","lon" )
Tps2 = right_join(Tps, ANVIO_cat2, by=c('Bin_Id'))
Tps2$lifestyle[Tps2 $lifestyle == "alga"] <- "algae"

Ctest  = cor.test(Tps2$Genome_size, Tps2 $rowSums.mat., method = "pearson")


hr1 = ggplot(Tps2 ,aes(x= Genome_size/1000000,y=rowSums.mat.))+geom_point(aes(color=group))+geom_smooth(method=lm,color='black',size=0.5)+ylab(paste('nr of ', COGselect_grep,sep=''))+ ggtitle(paste('pearsons cor ', round(Ctest$estimate,4),sep='')) +theme_classic() + theme(plot.title = element_text(size = 10, face = "bold"))

#hr2 = ggplot(Tps2 ,aes(x= Genome_size/1000000,y=rowSums.mat.))+geom_point(aes(color=lifestyle))+geom_smooth(method=lm,color='black',size=0.5)+ylab(paste('nr of ', COGselect_grep,sep=''))+theme_classic()
#hr3 = ggplot(Tps2 ,aes(x= Genome_size/1000000,y=rowSums.mat.,color=group))+geom_point()+geom_smooth(method=lm,fullrange=TRUE,se=FALSE,size=0.5)+ylab(paste('nr of ', COGselect_grep,sep=''))+theme_classic()
#hr4 = ggplot(Tps2 ,aes(x= Genome_size/1000000,y=rowSums.mat.,color=lifestyle))+geom_point()+geom_smooth(method=lm,fullrange=TRUE,se=FALSE,size=0.5)+ylab(paste('nr of ', COGselect_grep,sep=''))+theme_classic()

Tps2 $group <- factor(Tps2 $group, levels = c("adhaerens","algicola","lipo","psychro","antarcticus","hydrocarbo","excellens","vinifirmus","other"))

my_comparisons <- list( c("hydrocarbo", "lipo"), c("hydrocarbo", "algicola"), c("hydrocarbo", "antarcticus") )
p1 = ggplot(Tps2 ,aes(x= group,y=rowSums.mat./(Genome_size/1000000)))+geom_boxplot(
         )+xlab('lineage')+ylab(paste(COGselect_grep, ' related \n (genes per Mbp)',sep=''))+ theme_classic()+theme(axis.text.x=element_text(angle=45)) + stat_compare_means(comparisons = my_comparisons) + stat_compare_means()

my_comparisons <- list( c("Sediment", "Phototroph"),c("Photic", "Oil") )
p2 = ggplot(Tps2 ,aes(x= SS,y=rowSums.mat./(Genome_size/1000000)))+geom_boxplot(
         )+xlab('lineage')+ylab(paste(COGselect_grep, ' related \n (genes per Mbp)',sep=''))+ theme_classic()+theme(axis.text.x=element_text(angle=45)) + stat_compare_means(comparisons = my_comparisons) + stat_compare_means()

my_comparisons <- list( c("denitrifyer", "non-denitrifyer"),c("Photic", "Oil") )
p3 = ggplot(Tps2 ,aes(x= DNF2$DNF,y=rowSums.mat./(Genome_size/1000000)))+geom_boxplot(
         )+xlab('lineage')+ylab(paste(COGselect_grep, ' related \n (genes per Mbp)',sep=''))+ theme_classic()+theme(axis.text.x=element_text(angle=45)) + stat_compare_means(comparisons = my_comparisons) + stat_compare_means()


my_comparisons <- list( c("algae", "polar"),c("algae", "other"),c("polar",'other') )
p4 = ggplot(Tps2 ,aes(x= lifestyle,y=rowSums.mat./(Genome_size/1000000)))+geom_boxplot(
         )+xlab('lineage')+ylab(paste(COGselect_grep, ' related \n (genes per Mbp)',sep=''))+ theme_classic()+theme(axis.text.x=element_text(angle=45)) + stat_compare_means(comparisons = my_comparisons) + stat_compare_means()



#a <- aov(rowSums.mat./(Genome_size/1000000)~ group, data= Tps2)
#tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
#p1 = ggplot(Tps2 ,aes(x= group,y=rowSums.mat./(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab(paste(COGselect_grep, ' related \n (genes per Mbp)',sep=''))+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(Tps2$rowSums.mat./(Tps2$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

#a <- aov(rowSums.mat./(Genome_size/1000000)~ group, data= Tps2)
#tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
#p2 = ggplot(Tps2 ,aes(x= group,y=rowSums.mat.))+geom_boxplot()+xlab('lineage')+ylab(paste(COGselect_grep, ' related',sep=''))+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(Tps2$rowSums.mat.), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

multiplot(hhist1 ,hr1,p1,p2,p3,p4,cols=3)



mat_HEAT = mat[,1:ncol(mat)]
mat_HEAT2 = mat_HEAT
#mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]

mat_HEAT2_ordered = as.data.frame.matrix(mat_HEAT2_ordered)
heatmap.2(as.matrix(mat_HEAT2_ordered),trace="none"
,col = inferno(75),Rowv= rep_tree_d,srtCol=45,cexRow=0.5,cexCol=0.6,main= COGselect_grep)

distance.col = dist(t(as.matrix(mat_HEAT2_ordered)), method = "euclidean")
cluster.col = hclust(distance.col, method = "average")
mat_HEAT2_ordered[mat_HEAT2_ordered == 0] <- NaN

heat2 = heatmap.2( as.matrix(mat_HEAT2_ordered),
           col = inferno(75),Rowv= rep_tree_d,srtCol=45,
           trace = "none", 
           na.color="grey90",
           cexRow = 0.5, cexCol = 0.6,Colv=as.dendrogram(cluster.col),
           margins=c(12,12),keysize=0.75,main=COGselect_grep)





#ggplot(Tps2 ,aes(x= group,y=GC))+geom_boxplot()+xlab('lineage')+ylab('%GC')+theme_classic()+theme(axis.text.x=element_text(angle=45))
#ggplot(Tps2 ,aes(x= lifestyle,y=GC))+geom_boxplot()+xlab('lineage')+ylab('%GC')+theme_classic()+theme(axis.text.x=element_text(angle=45))


##########################
#KEGG
# KEGG and COG intervoncersion in based on the ko2kegg file on the KEGG ftp server
# file download the file from http://www.genome.jp/kegg/files/ko2cog.xl##########################
##########################
##########################################
##########################################

# Run this only if you need to, as it takes couple of hours
#extPath = c()
#for (i in ko2COG2$KO){
#	tryCatch({
#		print(i)
#       lkp <- keggGet(i)
# 	    infExtract=c(KO = i, KEGG_NAME = lkp[[1]]$NAME,KEGG_DEFINITION = lkp[[1]]$DEFINITION)
#        extPath = rbind(extPath, infExtract)
#     }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#}
# write.table(data.frame(extPath),file = 'KO_annotation2.txt', sep = "\t")
##########################################
##########################################


install.packages("KEGGREST",dependencies=TRUE)
library(KEGGREST)
library(stringr)

ko2COG = read.table("~/ko2cog.xl.txt",header=TRUE,sep='\t')
ko2COG2 = ko2COG
ko2COG2$COG_FUNCTION_ACC = gsub(" ","-",str_sub(ko2COG2$COG,6,-2))

ko2COG_annotated2 = read.table("~/KO_annotation2.txt",header=TRUE,sep='\t')
ko2COG_annotated3 = merge(ko2COG_annotated2, ko2COG2,by='KO')
ANVIOKEGG = merge(ANVIO, ko2COG_annotated3,by='COG_FUNCTION_ACC')



org <- keggList("pathway")
head(org)
org <- keggList("organism")
head(org)
org <- keggList("ko")
head(org)
queryables <- c(listDatabases(), org[,1], org[,2])
# the following selects a KEGG pathway and merges, and beautifies the associated COG identifiers into a single dataframe
query <- keggGet("ko00020") 
query <- keggGet("ko02010")

PathwKEGGset = data.frame(query[[1]]$ORTHOLOGY)

names(PathwKEGGset) = c("kegg_function")

PathwKEGGset$KO = rownames(PathwKEGGset) 
PathMerged = merge(PathwKEGGset, ko2COG,by='KO')

PathMerged$cog = gsub(" ","-",str_sub(PathMerged$COG,6,-2))

mat = t(table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% PathMerged$cog)[,c('COG_FUNCTION','genome_name')])))

mat_HEAT = mat[,1:ncol(mat)]
mat_HEAT2 = mat_HEAT
#mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]

mat_HEAT2_ordered = as.data.frame.matrix(mat_HEAT2_ordered)
heatmap.2(as.matrix(mat_HEAT2_ordered),trace="none"
,col = inferno(75),Rowv= rep_tree_d,srtCol=45,cexRow=0.5,cexCol=0.6,main='fructose')


test = c("K01535 ","K00329","K00330","K00331","K00332","K00333","K00334","K00335","K00336","K00337","K00338","K00339","K00340","K00341","K00342","K00343","K03878","K03879","K03880","K03881","K03882","K03883","K03884","K03934","K03935","K03936","K03940","K03941","K03942","K03943","K05572","K05573","K05574","K05575","K05576","K05577","K05578","K05579","K05580","K05581","K05582","K05583","K05584","K05585","K05586","K05587","K05588","K13378","K13380","K15863")
test = c('K13811','K00958','K00955','K00957','K00956','K13811','K00860 ','K00955','K00390','K00380','K00381','K00392','K02048','K02046','K02047','K02045')

PathMerged = merge(data.frame(KO = test), ko2COG,by='KO')
PathMerged$cog = gsub(" ","-",str_sub(PathMerged$COG,6,-2))
mat = t(table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% PathMerged$cog)[,c('COG_FUNCTION','genome_name')])))

mat_HEAT = mat[,1:ncol(mat)]
mat_HEAT2 = mat_HEAT
#mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]

mat_HEAT2_ordered = as.data.frame.matrix(mat_HEAT2_ordered)
heatmap.2(as.matrix(mat_HEAT2_ordered),trace="none"
,col = inferno(75),Rowv= rep_tree_d,srtCol=45,cexRow=0.5,cexCol=0.6,main='sulfur')


##############################################################################
##############################################################################


test = c('K13811','K00958','K00955','K00957','K00956','K13811','K00860 ','K00955','K00390','K00380','K00381','K00392','K02048','K02046','K02047','K02045')

mat = t(table(droplevels(subset(ANVIOKEGG, KO  %in% test)[,c('KEGG_DEFINITION','genome_name')])))

mat_HEAT = mat[,1:ncol(mat)]
mat_HEAT2 = mat_HEAT
#mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]

mat_HEAT2_ordered = as.data.frame.matrix(mat_HEAT2_ordered)
heatmap.2(as.matrix(mat_HEAT2_ordered),trace="none"
,col = inferno(75),Rowv= rep_tree_d,srtCol=45,cexRow=0.5,cexCol=0.6,main='sulfur')

##############################################################################
##############################################################################



KeggPathways=data.frame(keggList("pathway"))
names(KeggPathways) = c("pathway")
KeggPathways$pthw_accession = rownames(KeggPathways)

extPath = c()

for (i in KeggPathways$pthw_accession){
  lkp <- keggGet(i)
  infExtract=c(class = lkp[[1]]$CLASS,pathway_map = lkp[[1]]$PATHWAY_MAP,KO_pathway=lkp[[1]]$KO_PATHWAY)
  extPath = rbind(extPath, infExtract)
}

names(extPath)=c("class","Pathway_MAP","KO_PATHWAY")
KeggPathways = cbind(KeggPathways, extPath)


subset( KeggPathways, class %in% c("Metabolism; Energy metabolism","Cellular Processes; Cell motility","Environmental Information Processing; Membrane transport"))

#for (i in subset( KeggPathways, class %in% c("Metabolism; Carbohydrate metabolism"))$KO_pathway){
for (i in KeggPathways$KO_pathway){
tryCatch({
		print(i)
		query <- keggGet(i)
PathwKEGGset = data.frame(query[[1]]$ORTHOLOGY)
names(PathwKEGGset) = c("kegg_function")
PathwKEGGset$KO = rownames(PathwKEGGset) 
PathMerged = merge(PathwKEGGset, ko2COG,by='KO')
PathMerged$cog = gsub(" ","-",str_sub(PathMerged$COG,6,-2))

mat = t(table(droplevels(subset(ANVIOKEGG, KEGG_DEFINITION  %in% PathMerged$kegg_function)[,c('KEGG_DEFINITION','genome_name')])))

if(dim(mat)[1]==60){
mat_HEAT = mat[,1:ncol(mat)]
mat_HEAT2 = mat_HEAT
#mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]

mat_HEAT2_ordered = as.data.frame.matrix(mat_HEAT2_ordered)

distance.col = dist(t(as.matrix(mat_HEAT2_ordered)), method = "euclidean")
cluster.col = hclust(distance.col, method = "average")


mat_HEAT2_ordered[mat_HEAT2_ordered == 0] <- NaN


pdf(paste(i,".pdf", sep = ""))

heatmap.2( as.matrix(mat_HEAT2_ordered),
           col = inferno(75),Rowv= rep_tree_d,srtCol=45,
           trace = "none", 
           na.color="grey90",
           cexRow = 0.5, cexCol = 0.6,Colv=as.dendrogram(cluster.col),
           main=paste(i,query[[1]]$NAME,'\n',query[[1]]$CLASS,sep = " "),margins=c(12,12),keysize=0.75)
dev.off()

}}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}







##########################
# boxplots and significance on whatever
##########################

generate_label_df <- function(HSD, flev){
 # Extract labels and factor levels from Tukey post-hoc 
 Tukey.levels <- HSD[[flev]][,4]
 Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
 plot.labels <- names(Tukey.labels[['Letters']])

 # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
 # upper quantile and label placement
    boxplot.df <- ddply(d, flev, function (x) max(fivenum(x$y)) + 0.2)

 # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
     stringsAsFactors = FALSE)

 # Merge it with the labels
   labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)

return(labels.df)
}





COG_list = c("D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W" )
sum=0
mat = t(table(droplevels(ANVIO[grep("C", ANVIO $COG_CATEGORY),c('COG_FUNCTION','genome_name')])))
CpcogFull = data.frame(rowSums(mat))
for(i in COG_list){
	mat = t(table(droplevels(ANVIO[grep(i, ANVIO$COG_CATEGORY),c('COG_FUNCTION','genome_name')])))
	CpcogFull = cbind(CpcogFull,rowSums(mat))
}
names(CpcogFull)= c( "C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W" )
CpcogFull = cbind(CpcogFull, Bin_Id = rownames(CpcogFull))
CpcogFull = right_join(CpcogFull, CheckM2, by=c('Bin_Id'))
CpcogFull = right_join(CpcogFull, ANVIO_cat2, by=c('Bin_Id'))

CpcogFull $group <- factor(CpcogFull $group, levels = c("adhaerens","algicola","lipo","psychro","antarcticus","hydrocarbo","excellens","vinifirmus","other"))


a <- aov(y~lev, data=d)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
generate_label_df(tHSD, 'lev')

group = CpcogFull$group
y = CpcogFull$C/(CpcogFull$Genome_size/1000000)
a <- aov(C/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
 generate_label_df(tHSD, 'group')

ggplot(CpcogFull, aes(x= group, y=C/(Genome_size/1000000))) + geom_boxplot() +
  geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(y), label = labels))



max(CpcogFull$C/(CpcogFull$Genome_size/1000000))
max(CpcogFull$D/(CpcogFull$Genome_size/1000000))
max(CpcogFull$E/(CpcogFull$Genome_size/1000000))
max(CpcogFull$F/(CpcogFull$Genome_size/1000000))


a <- aov(C/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p1 = ggplot(CpcogFull ,aes(x= group,y=C/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('C \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$C/(CpcogFull$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(D/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p2 = ggplot(CpcogFull ,aes(x= group,y=D/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('D \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$D/(CpcogFull$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(E/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p3 = ggplot(CpcogFull ,aes(x= group,y=E/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('E \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$E/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(F/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p4 = ggplot(CpcogFull ,aes(x= group,y=F/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('F \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$F/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(G/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p5 = ggplot(CpcogFull ,aes(x= group,y=G/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('G \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$G/(CpcogFull$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(H/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p6 = ggplot(CpcogFull ,aes(x= group,y=H/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('H \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$H/(CpcogFull$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(I/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p7 = ggplot(CpcogFull ,aes(x= group,y=I/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('I \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$I/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(J/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p8 = ggplot(CpcogFull ,aes(x= group,y=J/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('J \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$J/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))


a <- aov(K/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p9 = ggplot(CpcogFull ,aes(x= group,y=K/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('K \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$K/(CpcogFull$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(L/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p10 = ggplot(CpcogFull ,aes(x= group,y=L/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('L \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$L/(CpcogFull$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(M/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p11 = ggplot(CpcogFull ,aes(x= group,y=M/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('M \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$M/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(N/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p12 = ggplot(CpcogFull ,aes(x= group,y=N/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('N \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$N/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(O/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p13 = ggplot(CpcogFull ,aes(x= group,y=O/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('O \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$O/(CpcogFull$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(P/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p14 = ggplot(CpcogFull ,aes(x= group,y=P/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('P \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$P/(CpcogFull$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(Q/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p15 = ggplot(CpcogFull ,aes(x= group,y=Q/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('Q \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$Q/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(R/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p16 = ggplot(CpcogFull ,aes(x= group,y=R/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('R \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$R/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(S/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p17 = ggplot(CpcogFull ,aes(x= group,y=O/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('O \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$O/(CpcogFull$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(T/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p18 = ggplot(CpcogFull ,aes(x= group,y=T/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('T \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$T/(CpcogFull$Genome_size/1000000)), label = labels))+ theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(U/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p19 = ggplot(CpcogFull ,aes(x= group,y=U/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('U \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$U/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(V/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p20 = ggplot(CpcogFull ,aes(x= group,y=V/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('V \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$V/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))

a <- aov(W/(Genome_size/1000000)~ group, data= CpcogFull)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
p21 = ggplot(CpcogFull ,aes(x= group,y=W/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('W \n (genes per Mbp)')+ geom_text(data = generate_label_df(tHSD, 'group'), aes(x = plot.labels, y = max(CpcogFull$W/(CpcogFull$Genome_size/1000000)), label = labels))+theme_classic()+theme(axis.text.x=element_text(angle=45))


multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,cols=6)




p5 = ggplot(CpcogFull ,aes(x= group,y=G/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('G \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))


p6 = ggplot(CpcogFull ,aes(x= group,y=H/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('H \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p7 = ggplot(CpcogFull ,aes(x= group,y=I/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('I \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p8 = ggplot(CpcogFull ,aes(x= group,y=J/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('J \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p9 = ggplot(CpcogFull ,aes(x= group,y=K/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('K \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p10 = ggplot(CpcogFull ,aes(x= group,y=L/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('L \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p11 = ggplot(CpcogFull ,aes(x= group,y=M/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('M \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p12 = ggplot(CpcogFull ,aes(x= group,y=N/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('N \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p13 = ggplot(CpcogFull ,aes(x= group,y=O/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('O \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p14 = ggplot(CpcogFull ,aes(x= group,y=P/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('P \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p15 = ggplot(CpcogFull ,aes(x= group,y=Q/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('Q \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p16 = ggplot(CpcogFull ,aes(x= group,y=R/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('R \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p17 = ggplot(CpcogFull ,aes(x= group,y=S/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('S \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p18 = ggplot(CpcogFull ,aes(x= group,y=T/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('T \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p19 = ggplot(CpcogFull ,aes(x= group,y=U/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('U \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p20 = ggplot(CpcogFull ,aes(x= group,y=V/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('V \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))
p21 = ggplot(CpcogFull ,aes(x= group,y=W/(Genome_size/1000000)))+geom_boxplot()+xlab('lineage')+ylab('W \n (genes per Mbp)')+theme_classic()+theme(axis.text.x=element_text(angle=45))

multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,cols=6)





plasmids = data.frame(rowSums(mat))
plasmids = data.frame(rowSums(mat), Bin_Id = rownames(plasmids))
Tps2 = right_join(Tps2, plasmids, by=c('Bin_Id'))

hh1 =ggplot(Tps2 ,aes(x= rowSums.mat..y,y=rowSums.mat..x))+geom_point(aes(color=group))+geom_smooth(method=lm,color='black',size=0.5)+ylab('nr of transposases')+ xlab('plasmid-related-genes') + theme_classic()
hh2= ggplot(Tps2 ,aes(x= rowSums.mat..y,y=rowSums.mat..x,color=lifestyle)) +geom_point(aes(color=lifestyle))+geom_smooth(method=lm,color='black',size=0.5)+ylab('nr of transposases')+xlab('plasmid-related-genes')+theme_classic()
multiplot(hh1, hh2,cols=2)
cor.test(Tps2$rowSums.mat..y, Tps2 $rowSums.mat..x, method = "pearson")

h21 =ggplot(Tps2 ,aes(x= GC,y=rowSums.mat..x))+geom_point(aes(color=group))+geom_smooth(method=lm,color='black',size=0.5)+ylab('nr of transposases')+ xlab('GC') + theme_classic()
h22= ggplot(Tps2 ,aes(x= GC,y=rowSums.mat..x,color=lifestyle)) +geom_point(aes(color=lifestyle))+geom_smooth(method=lm,color='black',size=0.5)+ylab('nr of transposases')+xlab('GC')+theme_classic()
multiplot(h21, h22,cols=2)

plot(Tps2$rowSums.mat., Tps2$Genome_size)

mat = t(table(droplevels(subset(ACCESmetadata2, COG_FUNCTION_ACC  %in% COGselect)[,c('COG_FUNCTION_ACC','genome_name')])))

mat = t(table(droplevels(subset(ACCESmetadata2, COG_FUNCTION  %in% c('DNA-binding'))[,c('COG_FUNCTION','genome_name')])))


mat_HEAT = mat[,2:ncol(mat)]
mat_HEAT2 = mat_HEAT
#mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]

mat_HEAT2_ordered = as.data.frame.matrix(mat_HEAT2_ordered)
annot_df <- data.frame(adhaerens = subset(STAMP1, COG  %in% colnames(mat_HEAT2_ordered))[,7],algicola = subset(STAMP2, COG  %in% colnames(mat_HEAT2_ordered))[,7],antarcticus = subset(STAMP3, COG  %in% colnames(mat_HEAT2_ordered))[,7],excellens = subset(STAMP4, COG  %in% colnames(mat_HEAT2_ordered))[,7],hydrocarb = subset(STAMP5, COG  %in% colnames(mat_HEAT2_ordered))[,7], psychro = subset(STAMP6, COG  %in% colnames(mat_HEAT2_ordered))[,7],vinifirmus = subset(STAMP7, COG  %in% colnames(mat_HEAT2_ordered))[,7])


annot.df = annot_df>0.05
annot.df[annot.df == TRUE] <- '#808080' 
annot.df[annot.df == FALSE] <- '#FFD700'

heatmap.plus(as.matrix(mat_HEAT2_ordered),trace="none",ColSideColors= annot.df
,col = inferno(75),Rowv= rep_tree_d,srtCol=45,cexRow=0.3,cexCol=0.3,main='highest contiruting')

heatmap.2(as.matrix(mat_HEAT2_ordered),trace="none",ColSideColors= annot.df[,1]
,col = inferno(max(mat_HEAT2_ordered)*2),Rowv= rep_tree_d,srtCol=45,cexRow=0.3,cexCol=0.3,main='highest contiruting')

heatmap.2(as.matrix(mat_HEAT2_ordered),trace="none"
,col = inferno(75),Rowv= rep_tree_d,srtCol=45,cexRow=0.3,cexCol=0.3,main='pilus')




mat = t(table(droplevels(subset(ACCESmetadata2, COG_FUNCTION_ACC  %in% c(as.character(STAMP1[STAMP1[,7]<0.01,]$COG),as.character(STAMP2[STAMP2[,7]<0.01,]$COG),as.character(STAMP3[STAMP3[,7]<0.01,]$COG),as.character(STAMP4[STAMP4[,7]<0.01,]$COG),as.character(STAMP5[STAMP5[,7]<0.01,]$COG),as.character(STAMP6[STAMP6[,7]<0.01,]$COG),as.character(STAMP7[STAMP7[,7]<0.01,]$COG),as.character(STAMP8[STAMP8[,7]<0.01,]$COG)
))[,c('combined','genome_name')])))  



mat = t(table(droplevels(subset(ACCESmetadata2, COG_FUNCTION_ACC  %in% c(as.character(STAMP1[STAMP1[,7]<0.01,]$COG),as.character(STAMP2[STAMP2[,7]<0.01,]$COG),as.character(STAMP3[STAMP3[,7]<0.01,]$COG),as.character(STAMP4[STAMP4[,7]<0.01,]$COG),as.character(STAMP5[STAMP5[,7]<0.01,]$COG),as.character(STAMP6[STAMP6[,7]<0.01,]$COG),as.character(STAMP7[STAMP7[,7]<0.01,]$COG),as.character(STAMP8[STAMP8[,7]<0.01,]$.)
))[,c('COG_FUNCTION_ACC','genome_name')])))

mat_HEAT = mat[,2:ncol(mat)]
mat_HEAT2 = mat_HEAT
#mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]

mat_HEAT2_ordered = as.data.frame.matrix(mat_HEAT2_ordered)
annot_df <- data.frame(adhaerens = subset(STAMP1, COG  %in% colnames(mat_HEAT2_ordered))[,7],algicola = subset(STAMP2, COG  %in% colnames(mat_HEAT2_ordered))[,7],antarcticus = subset(STAMP3, COG  %in% colnames(mat_HEAT2_ordered))[,7],excellens = subset(STAMP4, COG  %in% colnames(mat_HEAT2_ordered))[,7],lipo = subset(STAMP8, .  %in% colnames(mat_HEAT2_ordered))[,7],hydrocarb = subset(STAMP5, COG  %in% colnames(mat_HEAT2_ordered))[,7], psychro = subset(STAMP6, COG  %in% colnames(mat_HEAT2_ordered))[,7],vinifirmus = subset(STAMP7, COG  %in% colnames(mat_HEAT2_ordered))[,7])


annot.df = annot_df>0.05
annot.df[annot.df == TRUE] <- '#808080' 
annot.df[annot.df == FALSE] <- '#FFD700'

heatmap.plus(as.matrix(mat_HEAT2_ordered),trace="none",ColSideColors= annot.df
,col = inferno(75),Rowv= rep_tree_d,srtCol=45,cexRow=0.3,cexCol=0.3,main='highest contiruting')

heatmap.2(as.matrix(mat_HEAT2_ordered),trace="none"
,col = inferno(75),Rowv= rep_tree_d,srtCol=45,cexRow=0.3,cexCol=0.3,main='highest contiruting')

annot_df <- data.frame(adhaerens = STAMP1[,7],algicola = STAMP2[,7],antarcticus = STAMP3[,7],excellens = STAMP4[,7],hydrocarb = STAMP5[,7], psychro = STAMP6[,7],vinifirmus = STAMP7[,7])



mat = t(table(droplevels(subset(ACCESmetadata2, COG_FUNCTION_ACC  %in% c(as.character(STAMP1[STAMP1[,7]<0.00000000001,]$COG),as.character(STAMP2[STAMP2[,7]<0.00000000001,]$COG),as.character(STAMP3[STAMP3[,7]<0.00000000001,]$COG),as.character(STAMP4[STAMP4[,7]<0.00000000001,]$COG),as.character(STAMP5[STAMP5[,7]<0.00000000001,]$COG),as.character(STAMP6[STAMP6[,7]<0.00000000001,]$COG),as.character(STAMP7[STAMP7[,7]<0.00000000001,]$COG)
))[,c('COG_FUNCTION_ACC','genome_name')])))

mat = t(table(droplevels(subset(ACCESmetadata2, COG_FUNCTION_ACC  %in% c(as.character(STAMP1[STAMP1[,7]<0.0000000000000001,]$COG),as.character(STAMP2[STAMP2[,7]<0.0000000000000001,]$COG),as.character(STAMP3[STAMP3[,7]<0.0000000000000001,]$COG),as.character(STAMP4[STAMP4[,7]<0.0000000000000001,]$COG),as.character(STAMP5[STAMP5[,7]<0.0000000000000001,]$COG),as.character(STAMP6[STAMP6[,7]<0.0000000000000001,]$COG),as.character(STAMP7[STAMP7[,7]<0.0000000000000001,]$COG)
))[,c('combined','genome_name')])))


mat = t(table(droplevels(subset(ACCESmetadata2, COG_FUNCTION_ACC  %in% c(as.character(STAMP1[STAMP1[,7]<0.01,]$COG),as.character(STAMP2[STAMP2[,7]<0.01,]$COG),as.character(STAMP3[STAMP3[,7]<0.01,]$COG),as.character(STAMP4[STAMP4[,7]<0.01,]$COG),as.character(STAMP5[STAMP5[,7]<0.01,]$COG),as.character(STAMP6[STAMP6[,7]<0.01,]$COG),as.character(STAMP7[STAMP7[,7]<0.01,]$COG)
))[,c('combined','genome_name')])))      

 

test <- read.table("~/DATA/SampleAnnotatedHeatmapDataNorm.txt",row.names = 1)
# control <- test[grep("Control",row.names(test)),]
# treatment_3weeks <- test[grep("Treatment_3weeks",row.names(test)),]
# GSM_data <- test[grep("GSM[1][0-2]",row.names(test)),]





STAMP1[STAMP1[3]>STAMP1[3]]


DEDUPLICATED =  subset(ANVIO, !duplicated(COG_FUNCTION_ACC))
colnames(DEDUPLICATED)=c("unique_id","gr","bin_name","genome_name","gene_callers_id","COG_CATEGORY_ACC","COG_CATEGORY","COG","COG_FUNCTION","aa_sequence" )

JOINED <- right_join(DEDUPLICATED, STAMP2[STAMP2[,2]>STAMP2[,4] & STAMP2[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/algicola_enriched_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP2[STAMP2[,2]<STAMP2[,4] & STAMP2[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/algicola_depleted_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP1[STAMP1[,2]>STAMP1[,4] & STAMP1[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/adhaerens_enriched_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP1[STAMP1[,2]<STAMP1[,4] & STAMP1[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/adhaerens_depleted_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP3[STAMP3[,2]>STAMP3[,4] & STAMP3[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/antarcticus_enriched_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP3[STAMP3[,2]<STAMP3[,4] & STAMP3[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/antarcticus_depleted_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP4[STAMP4[,2]>STAMP4[,4] & STAMP4[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/excellens_enriched_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP4[STAMP4[,2]<STAMP4[,4] & STAMP4[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/excellens_depleted_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP5[STAMP5[,2]>STAMP5[,4] & STAMP5[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/hydrocarb_enriched_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP5[STAMP5[,2]<STAMP5[,4] & STAMP5[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/hydrocarb_depleted_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP6[STAMP6[,2]>STAMP6[,4] & STAMP6[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/psychro_enriched_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP6[STAMP6[,2]<STAMP6[,4] & STAMP6[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/psychro_depleted_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP7[STAMP7[,2]>STAMP7[,4] & STAMP7[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/vinifirmus_enriched_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP7[STAMP7[,2]<STAMP7[,4] & STAMP7[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/vinifirmus_depleted_0.01.txt",sep="\t")

JOINED <- right_join(DEDUPLICATED, STAMP8[STAMP8[,2]>STAMP8[,4] & STAMP8[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/lipo_enriched_0.01.txt",sep="\t")
JOINED <- right_join(DEDUPLICATED, STAMP8[STAMP8[,2]<STAMP8[,4] & STAMP8[,7]<0.01,], by=c("COG"))
write.table(JOINED,"~/lipo_depleted_0.01.txt",sep="\t")


# a function to assign colors based on treatment time 
# http://www.rapidtables.com/web/color/RGB_Color.htm
condition_colors <- unlist(lapply(rownames(test),function(x){
  if(grepl("Treatment",x)) '#FFC0CB' #pink
  else if(grepl('Control',x)) '#808080' #grey
  
}))

# I like to have a line just to assign input in one step
input <- as.matrix(t(test))

heatmap.2(input, trace="none", density="none", col=bluered(20), cexRow=1, cexCol=0.2, margins = c(20,13),
          ColSideColors=condition_colors, scale="row")

treatment_times <- c(0,1,3,8,24)
treatment_color_options <- c(brewer.pal(5, "Set1"))

treatment_colors <- unlist(lapply(rownames(test),function(x){
  for(t_time in 1:length(treatment_times)){
    if(grepl(paste("_",treatment_times[t_time],"weeks",sep=""),x)) return(treatment_color_options[t_time])
  }
}))

myCols <- cbind(condition_colors,treatment_colors)
colnames(myCols)[1] <- "Condition"
colnames(myCols)[2] <- "Treatment Time"

# for exporting pdf or jpeg etc:
pdf(file="~/YourPath/AnnotatedHeatmap.pdf")  
par(cex.main=0.8,mar=c(1,1,1,1))
heatmap.plus(input, col=bluered(20),cexRow=1,cexCol=0.2, margins = c(20,13), main="Your Title",
ColSideColors=myCols)








col = list(adhaerens =  c("TRUE" = "red", "FALSE" = "white"),algicola = c("TRUE" = "red", "FALSE" = "white"), antarcticus = c("TRUE" = "red", "FALSE" = "white") )
deduplALLanvio = distinct(ANVIO, COG_CATEGORY_ACC, protein_cluster_id)



COG_list = c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W" )
sum=0
COGALL= c() 
COGCORE = c()
COGSINGL = c()
COGACCES = c()
for(i in COG_list){
	#nrCOGS = c(nrCOGS,nrow(ANVIOData[grep(i, ANVIOData$COG_CATEGORY_ACC),]))
	print(paste(i,nrow(deduplCOREanvio[grep(i, deduplCOREanvio $COG_CATEGORY_ACC),]),sep=' - '))
	sum = sum+nrow(deduplCOREanvio[grep(i, deduplCOREanvio $COG_CATEGORY_ACC),])
	COGCORE = c(COGCORE,nrow(deduplCOREanvio[grep(i, deduplCOREanvio $COG_CATEGORY_ACC),]))
	COGSINGL = c(COGSINGL,nrow(deduplSINGLanvio[grep(i, deduplSINGLanvio$COG_CATEGORY_ACC),]))
	COGACCES = c(COGACCES,nrow(deduplACCESanvio[grep(i, deduplACCESanvio $COG_CATEGORY_ACC),]))
	COGALL = c(COGALL,nrow(deduplALLanvio[grep(i, deduplALLanvio $COG_CATEGORY_ACC),]))
}

COGPAN = data.frame(COG = COG_list,CORE= COGCORE,ACCES= COGACCES , UNIQUE= COGSINGL,ALL= COGALL)

print(paste("conserved information strorage and processing: ",sum(subset(COGPAN,COG %in% c('J','K','L'))$CORE)/sum(subset(COGPAN,COG %in% c('J','K','L'))$ALL)))
print(paste("conserved cellular process and signaling: ",sum(subset(COGPAN,COG %in% c('D','V','W','T','M','N','U','O'))$CORE)/sum(subset(COGPAN,COG %in% c('D','V','W','T','M','N','U','O'))$ALL)))
print(paste("conserved metabolism: ",sum(subset(COGPAN,COG %in% c('C','G','E','F','H','I','P','Q'))$CORE)/sum(subset(COGPAN,COG %in% c('C','G','E','F','H','I','P','Q'))$ALL)))
print(paste("conserved poorly characterized: ",sum(subset(COGPAN,COG %in% c('R','S','-'))$CORE)/sum(subset(COGPAN,COG %in% c('R','S','-'))$ALL)))


COGPAN2 = data.frame(COG = COG_list,CORE= (COGCORE*100)/COGALL,ACCES= (COGACCES*100)/COGALL , UNIQUE= (COGSINGL*100)/COGALL)
levels(COGPAN2$COG) = c(levels(COGPAN2$COG),"-")

COGPAN2 = rbind(COGPAN2,c("-",nrow(deduplCOREanvio[deduplCOREanvio$COG_CATEGORY_ACC == '',])*100/nrow(deduplALLanvio[deduplALLanvio $COG_CATEGORY_ACC == '',]),nrow(deduplSINGLanvio[deduplSINGLanvio $COG_CATEGORY_ACC == '',])*100/nrow(deduplALLanvio[deduplALLanvio $COG_CATEGORY_ACC == '',]),nrow(deduplACCESanvio[deduplACCESanvio $COG_CATEGORY_ACC == '',])*100/nrow(deduplALLanvio[deduplALLanvio $COG_CATEGORY_ACC == '',])))

CPN = melt(COGPAN2,id=c('COG'))

CPN$COG2 <- factor(CPN$COG, rev(c("J", "K", "L", "D","V","W", "T", "M", "N","U", "O", "C", "G","E", "F", "H", "I","P", "Q", "R", "S","-")))

 ggplot(CPN, aes(x = COG2, y = as.numeric(as.character(value)), fill = factor(variable,levels=rev(c('CORE','ACCES','UNIQUE'))))) + 
  geom_bar(stat = "identity",width=0.8) + xlab("") + ylab("") + coord_flip() + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())

# Make figure showing proportions of COG assigned vs not assigned PCs per PAN genome class

NOCOG_proportions = c(nrow(deduplCOREanvio[deduplCOREanvio$COG_CATEGORY_ACC == '',])/ncol(COREAnv),nrow(deduplSINGLanvio[deduplSINGLanvio $COG_CATEGORY_ACC == '',])/ncol(ACCESLAnv),nrow(deduplSINGLanvio[deduplSINGLanvio $COG_CATEGORY_ACC == '',])/ncol(SINGLAnv))

Assignment = data.frame(PAN = c('CORE','ACCESSORY','UNIQUE'),NotAssigned= NOCOG_proportions, Assigned = c(1,1,1)-NOCOG_proportions)

ggplot(melt(Assignment), aes(x = factor(PAN,levels=rev(c('CORE','ACCESSORY','UNIQUE'))), y = as.numeric(as.character(value)), fill = factor(variable,levels=rev(c('Assigned','NotAssigned'))),label = round(as.numeric(as.character(value))*100,2))) + 
  geom_bar(stat = "identity",width=0.8)+geom_text(size = 3, position = position_stack(vjust = 0.5)) + xlab("") + ylab("") + coord_flip() + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())


###########
# obtain group specific core genomes
############



###########
# obtain group specific core genomes
############


for(phylogroup in unique(ANVIO_cat$group)){
	selection = ANVIO_cat[grep(phylogroup, ANVIO_cat$group),]
	metaSelection = AnvioData[rownames(AnvioData) %in% selection$name,]
	print(paste(phylogroup , nrow(metaSelection)))
	COREAnv = metaSelection[,colSums(metaSelection)==nrow(metaSelection) ]
	print(ncol(COREAnv))
	COREmetadata = ANVIO[ANVIO $protein_cluster_id %in% colnames(COREAnv),]
	COREmetadata2 = COREmetadata[COREmetadata$genome_name ==rownames(COREAnv)[grep(phylogroup,rownames(COREAnv))],]
	
	D = do.call(rbind, lapply(seq(nrow(COREmetadata2[, c('protein_cluster_id','aa_sequence')])), function(i) t(COREmetadata2[, c('protein_cluster_id','aa_sequence')][i, ])))

	write.table(D, file = paste('~/', phylogroup , '_CORE.fasta',sep=''),row.names = FALSE, col.names = FALSE, quote = FALSE)
	}
	
	
	
AnvioData[ANVIO$bin_name==i,]
ANVIO_cat[grep(i, ANVIO_cat$group),]
#----------------------------------------------
#	PAN GENOME CURVE
#----------------------------------------------
#---- species accumulation curve --------------

sp <- specaccum(AnvioData, 'random', permutations=100)
Pfam_sp <- specaccum(PfamData, 'random', permutations=100)

binomix2 <- binomixEstimate(AnvioData, K.range=2:15)
binomix3 <- binomixEstimate(PfamData, K.range=2:15)

summary(binomix2)
summary(binomix3)

summary(sp)
plot(sp, ci.type='poly', col='darkred', lwd=2, ci.lty=0, ci.col='darkred', xlab='Genomes', ylab='Proteins', main='Gene accumulation plot')
boxplot(sp, col='white', add=TRUE,cex=0.5,pch = 21) 
boxplot(Pfam_sp, col='darkgrey', add=TRUE,cex=0.5,pch = 21) 

mods <- fitspecaccum(sp, "arrh")
plot(mods, col="grey",,xlab='Genomes', ylab='Protein Clusters', main='Gene accumulation plot')
boxplot(sp, col = "white", border = "black", lty=1, cex=0.3, add= TRUE)
## Use nls() methods to the list of models
sapply(mods$models, AIC)

plot(sp, ci.type='poly', col='black', lwd=2, ci.lty=0, ci.col='lightskyblue3', xlab='Genomes', ylab='Proteins', main='accumulation plot')
boxplot(sp, col = "white", border = "darkblue", lty=1, cex=0.3, add= TRUE)
par(new = T)
plot(Pfam_sp, ci.type='poly', col='black', lwd=2, ci.lty=0, ci.col='pink', axes=F, xlab=NA, ylab=NA)
boxplot(Pfam_sp, col = "white", border = "darkred", lty=1, cex=0.3, add= TRUE)
axis(side = 4)
mtext(side = 4, line = 3, 'PfamDomains')


#---- Rarefraction curve --------------
#	Alternative way to plot Pangenome curve

rarefrac <- rarefaction(AnvioData,n.perm=100)
    plot(rarefrac)
 summary(rarefrac)

rarefrac_Pfam <- rarefaction(PfamData,n.perm=100)
    plot(rarefrac_Pfam)
 summary(rarefrac_Pfam)

#---- sopen or closed panGENOME --------------
# Estimating if the pan-genome is open or closed based on a Heaps law model
#This function is based on a Heaps law approach suggested by Tettelin et al (2008). The Heaps law model is fitted to the number of new gene clusters observed when genomes are ordered in a random way. The model has two parameters, an intercept and a decay parameter called alpha. If alpha>1.0 the pan-genome is closed, if alpha<1.0< span=""> it is open.
#----------------------------------------------

heaps(AnvioData,n.perm=1000)
heaps(PfamData,n.perm=1000)


#---- Pan-genome size --------------
#Chao - computes the Choa lower bound estimated number of gene clusters in a pan-genome

chao.pansize <- chao(AnvioData)
chao.Pfampansize <- chao(PfamData)

chao.pansize
chao.Pfampansize


#---- DISTANCES AND WEIGHTED DISTANCES --------------
#Manhattan and/or Jaccard distances between pan-genome profiles. Jaccard distance is based on similarity of shared genes, while Manhattan distance also includes similarity of lacking a certain gene (which is often recommended).(Snipen and Ussery, 2010)
#----------------------------------------------

#Jaccard distances based on panmatrix
Jdist.blast <- distJaccard(AnvioData)

#manhatten distamces 
Mdist.blast <- distManhattan(AnvioData)

#fluidityy
fluid.blast <- fluidity(AnvioData)


#----------------------------------------------
#This weighting function computes weights for gene cluster according to their distribution in a pan-genome. 
#When computing distances between genomes or a PCA, it is possible to give weights to the different gene clusters, emphasizing certain aspects. 
#As proposed by Snipen & Ussery (2010), we have implemented two types of weighting: The default ‘"shell"’ type means gene families occuring frequently in the genomes, denoted shell-genes, are given large weight (close to 1) while those occurring rarely are given small weight (close to 0). The opposite is the ‘"cloud"’ type of weighting. Genes observed in a minority of the genomes are referred to as cloud-genes. Presumeably, the ‘"shell"’ weighting will give distances/PCA reflecting a more long-term evolution, since emphasis is put on genes who have just barely diverged away from the core. The ‘"cloud"’ weighting emphasizes those gene clusters seen rarely. Genomes with similar patterns among these genes may have common recent history. A ‘"cloud"’ weighting typically gives a more erratic or ‘noisy’ picture than the ‘"shell"’ weighting. 
#----------------------------------------------

w <- geneWeights(AnvioData,type="shell")
w2 <- geneWeights(AnvioData,type="cloud")

Mdist.blast <- distManhattan(AnvioData,weights=w)
Jdist.blast <- distJaccard(AnvioData,weights = w)


m2 <- melt(as.matrix(Mdist.blast))[melt(upper.tri(as.matrix(Mdist.blast)))$value,]
m3 <- melt(as.matrix(Jdist.blast))[melt(upper.tri(as.matrix(Jdist.blast)))$value,]

constructedMatrix = data.frame(cbind(m2$value,m3$value))
colnames(constructedMatrix) = c("Manthattan","Jaccard")

ggplot(constructedMatrix, aes(Manthattan, Jaccard)) + geom_point() + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) 


# # Based on these results, we can reject the null hypothesis that these two matrices, spatial distance and ozone distance, are unrelated with alpha = .05. The observed correlation, r = 0.1636308, suggests that the matrix entries are positively associated.  So smaller differences in ozone are generally seen among pairs of stations that are close to each other than far from each other. Note that since this test is based on random permutations, the same code will always arrive at the same observed correlation but rarely the same p-value.
mantel.rtest(Mdist.blast, Jdist.blast, nrepet = 999)



datadist <- distManhattan(AnvioData)
shelldist <- distManhattan(AnvioData,weights=w)
clouddist <- distManhattan(AnvioData,weights=w2)

m1 <- melt(as.matrix(datadist))[melt(upper.tri(as.matrix(datadist)))$value,]
m2 <- melt(as.matrix(shelldist))[melt(upper.tri(as.matrix(shelldist)))$value,]
m3 <- melt(as.matrix(clouddist))[melt(upper.tri(as.matrix(clouddist)))$value,]

constructedMatrix = data.frame(cbind(m2$value,m3$value))
colnames(constructedMatrix) = c("shell","cloud")

ggplot(constructedMatrix, aes(shell, cloud)) + geom_point() + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) 

constructedMatrix = data.frame(cbind(m1$value,m2$value))
colnames(constructedMatrix) = c("pan","shell")

ggplot(constructedMatrix, aes(pan, shell)) + geom_point() + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T)  + geom_smooth(method="lm", formula= (y ~ exp(x)), se=FALSE, linetype = 1)

ggplot(constructedMatrix, aes(pan, shell)) + geom_point() + geom_smooth(method="lm", formula= (y ~ exp(x)), se=FALSE, linetype = 1)


mantel.rtest(shelldist, clouddist, nrepet = 999)
mantel.rtest(datadist, shelldist, nrepet = 999)


v <- geneWeights(PfamData,type="shell")
v2 <- geneWeights(PfamData,type="cloud")

PfMdist.blast <- distManhattan(PfamData,weights=v)
PfJdist.blast <- distJaccard(PfamData,weights = v)


#---- NETWORK TREE BASED ON WEIGHTED DISTANCES --------------
#FOR ANVIO GENE CONTENT
nnet <- neighborNet(Mdist.blast)
#par("mar" = rep(1, 4))
plot(nnet, "2D")
plot(nnet)

#FOR Pfam Domain content
Pnnet <- neighborNet(PfMdist.blast)
#par("mar" = rep(1, 4))
plot(Pnnet, "2D")
plot(Pnnet)


#---- PRINCIPAL COMPONENT ANALYSIS --------------
#	PCA analysis weighted or not
#----------------------------------------------
pca = panpca(AnvioData, scale = 0, weights = w)
pca2 = panpca(AnvioData, scale = 0, weights = w2)
pca3 = panpca(AnvioData, scale = 0)
pca4 = panpca(AnvioData, scale = 0)
pca4 = panpca(AnvioData, scale = 0)


pfpca = panpca(PfamData, scale = 0, weights = v)
pfpca2 = panpca(PfamData, scale = 0, weights = v2)

ANVIO_cat <- read.table('~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/ANVIO_CAT.txt',header=TRUE)

cols <- RColorBrewer::brewer.pal(length(unique(ANVIO_cat$group)), name = "RdYlBu")
#cols = c('black','red',"grey",'pink','yellow','cyan3','green','purple','orange')
ANVIO_cat$color <- factor(ANVIO_cat$group, labels = cols)

pcacolors <- ANVIO_cat$color
names(pcacolors) <- ANVIO_cat$name
names(pcacolors) <- rownames(AnvioData)

#labels <- ANVIO_cat$name
#names(labels) <- ANVIO_cat$name

plotScores(pca, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
plotScores(pca2, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
plotScores(pca3, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)

legend("topright", col=unique(ANVIO_cat$color), legend = unique(ANVIO_cat$group),
    pch = 20, bty='n', cex=.75)

plotScores(pca,col= pcacolors)

    
'COG0069-COG0070','COG0489-COG3536','COG0503'

layout3d(matrix(1:4, 2, 2),sharedMouse = TRUE)

plot3d(pca$Scores[,1:3], col= pcacolors, size=20, type='p')
grid3d(c("x", "y+", "z"))
#text3d(pca$Scores[,1], pca$Scores[,2], pca$Scores[,3],texts=c(rownames(pca$Scores)), cex= 0.7, pos=3)      
next3d()

plot3d(pca2$Scores[,1:3], col= pcacolors, size=20, type='p',main='test')
grid3d(c("x", "y+", "z"))
#text3d(pca2$Scores[,1], pca2$Scores[,2], pca2$Scores[,3],texts=c(rownames(pca2$Scores)), cex= 0.7, pos=3)
next3d()


plot3d(pfpca$Scores[,1:3], col= pcacolors, size=20, type='p')
grid3d(c("x", "y+", "z"))
#text3d(pfpca$Scores[,1], pfpca $Scores[,2], pfpca$Scores[,3],texts=c(rownames(pfpca$Scores)), cex= 0.7, pos=3) 


#looking at DUF group only
pfpca = panpca(PfamData[,grep("DUF", colnames(PfamData))], scale = 0, weights = v)
plot3d(pfpca$Scores[,1:3], col= pcacolors, size=20, type='p')
grid3d(c("x", "y+", "z"))
plotScores(pfpca, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
#looking at DUF group only

pfpca = panpca(PfamData[,grep("Phage", colnames(PfamData))], scale = 0, weights = v)
plot3d(pfpca$Scores[,1:3], col= pcacolors, size=20, type='p')
grid3d(c("x", "y+", "z"))

next3d()
plot3d(pfpca2$Scores[,1:3], col= pcacolors, size=20, type='p')
grid3d(c("x", "y+", "z"))
#text3d(pfpca2$Scores[,1], pfpca2 $Scores[,2], pfpca2$Scores[,3],texts=c(rownames(pfpca2$Scores)), cex= 0.7, pos=3)

   
rgl.postscript("PCA_deluxe_ran.pdf","pdf") 


COG_profile = ANVIO_cat

COG_profile$C = as.numeric(table(SLCTANV$genome))
p2 = ggplot(COG_profile, aes(x = name, y = C, fill= group)) + geom_boxplot() + xlab("pan-genome devision") + ylab("number of PCs") + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) + geom_hline(yintercept = unique(Core_count),linetype = 2)


par(mfrow=c(5,5), mar = c(5,5,5,5))
SLCTANV = ANVIO[grep("C", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16,main="C")
COG_profile$C = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("D", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16,xlim=c(-6,-6))
COG_profile$D = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("E", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$E = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("F", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$F = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("G", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$G = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("H", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$H = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("I", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$I = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("J", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$J = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("K", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$K = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("L", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$L = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("M", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$M = as.numeric(table(SLCTANV$genome))


SLCTANV = ANVIO[grep("N", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$N = as.numeric(table(SLCTANV$genome))


SLCTANV = ANVIO[grep("O", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$O = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("P", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$P = as.numeric(table(SLCTANV$genome))

SLCTANV = ANVIO[grep("Q", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$Q = as.numeric(table(SLCTANV$genome))


SLCTANV = ANVIO[grep("R", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$R = as.numeric(table(SLCTANV$genome))


SLCTANV = ANVIO[grep("S", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$S = as.numeric(table(SLCTANV$genome))


SLCTANV = ANVIO[grep("T", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$T = as.numeric(table(SLCTANV$genome))


SLCTANV = ANVIO[grep("U", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$U = as.numeric(table(SLCTANV$genome))


SLCTANV = ANVIO[grep("V", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$V = as.numeric(table(SLCTANV$genome))


SLCTANV = ANVIO[grep("W", ANVIO$COG_CATEGORY), ]
pca4 = panpca(AnvioData[,unique(SLCTANV$protein_cluster_id)], scale = 0)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
COG_profile$W = as.numeric(table(SLCTANV$genome))

subselection = data.frame(cbind(as.character(SLCTANV$COG_FUNCTION),as.character(SLCTANV$genome_name)))
colnames(subselection)=c('COG_FUNCTION','name')

dunno = data.frame((table(subselection)))
dunno2 = (merge(dunno, ANVIO_cat, by = 'name'))
ggplot(dunno2, aes(x= COG_FUNCTION,y=Freq)) + 
    geom_col(aes(fill=group)) +   coord_flip() + theme(text = element_text(size=5)) + scale_fill_manual(values = getPalette(9))

ggplot(subset(dunno2, COG_FUNCTION == 'Beta-lactamase_class_A'), aes(x= name,y=Freq)) + 
    geom_col() +   coord_flip() + theme(text = element_text(size=10)) + scale_x_discrete(limits=tree2$tip.label[ordered_tips])


clade_order <- order.dendrogram(rep_tree_d)
clade_name <- labels(rep_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(AccesoryDF))

#SAVES AS PDF, NO VISUALISATION
COG_list = c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W" )
for(i in COG_list){
SLCTANV = ANVIO[grep(i, ANVIO$COG_CATEGORY), ]
subselection = data.frame(cbind(as.character(SLCTANV$COG_FUNCTION),as.character(SLCTANV$genome_name)))
colnames(subselection)=c('COG_FUNCTION','name')
dunno = data.frame((table(subselection)))
COG_HEAT = reshape(dunno, idvar = "name", timevar = "COG_FUNCTION", direction = "wide")
rownames(COG_HEAT) =  COG_HEAT$name
COG_HEAT = COG_HEAT[,2:ncol(COG_HEAT)]
COG_HEAT2 = COG_HEAT
COG_HEAT2[COG_HEAT2>1] <-1
COG_HEAT2_ordered <- COG_HEAT2[new_order,]

pdf(paste(i,"_heat.pdf",sep=''))
heatmap.2(as.matrix(t(COG_HEAT2_ordered)),trace="none",col =c('white',"black"),Colv= rep_tree_d,margins = c(1000/ncol(COG_HEAT2), 10),dendrogram="col",cexCol=0.5,cexRow=10/ncol(COG_HEAT2)*3,key=FALSE)
#dev.copy(pdf,paste(i,"_heat.pdf",sep=''))
dev.off()
}
COG_list = c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W" )
cFasta= c()

COG_list = c("J" )
for(i in COG_list){
SLCTANV = ANVIO[grep(i, ANVIO$COG_CATEGORY), ]
for(o in unique(SLCTANV[, "protein_cluster_id"])){
	SLCTANV = ANVIO[grep(o, ANVIO$protein_cluster_id), ]
	if( length(unique(SLCTANV$genome_name)) == 60 ){
		cFasta= rbind(cFasta,data.frame(SLCTANV$genome_name,SLCTANV$aa_sequence, SLCTANV$COG_FUNCTION))

}
}
}
write.csv(cFasta,"J_CORE.csv")
 
 
 
install.packages("seqRFLP")
library("seqRFLP")

df.fasta = dataframe2fas(data.frame(cFasta[1:60,c("SLCTANV.genome_name","SLCTANV.aa_sequence")]), file="df.fasta")
 aln2 <- muscle(aln)
aln = read.alignment('df.fasta',format = 'fasta')
aln = read.alignment('df.fasta',format = 'fasta')


par(mfrow=c(5,5), mar = c(5,5,5,5))
l <- list();
o <- 1

COG_list = c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W" )
for(i in COG_list){
  plt <- ggplot(COG_profile, aes_string(x= COG_profile$group, y = i)) + geom_boxplot()  + theme_classic() + theme(axis.title.x=element_blank())
  print(plt)
  l[[o]] <- plt
  o = o +1
}
layout <-  matrix(c(1:25), nrow = 5, byrow = TRUE)

multiplot(plotlist = l,layout=layout)


getPalette = colorRampPalette(brewer.pal(9,'Set1'))
COG_prof_long <- melt(COG_profile, id=c("name","group","lifestyle","color" )) 
g = ggplot(COG_prof_long,aes(x=name,y=value)) + geom_col(aes(fill=variable)) + coord_flip() + scale_fill_manual(values = getPalette(21))


RAxMLANVIO = read.tree("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/ANVIO_RAxML/RAxML_bipartitions.autosubst100")
RAxMLANVIORooted =  midpoint.root(RAxMLANVIO)
tree2 <- ladderize(RAxMLANVIORooted, right = FALSE)
RAxMLANVIORooted =  midpoint.root(RAxMLANVIO)
tree2 <- ladderize(RAxMLANVIORooted, right = FALSE)
tr = plot(tree2)
is_tip <- tree2$edge[,2] <= length(tree2$tip.label)
ordered_tips <- tree2$edge[is_tip, 2]
tree2$tip.label[ordered_tips]

PROF = ggplot(COG_prof_long, aes(x=name,y=value)) + 
    geom_col(aes(fill=variable)) + 
    xlab("Genome") +
    ylab("COG profile") + scale_fill_manual(values = getPalette(21)) +
    theme_classic() + scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()
PROF

COG_profile_perc = cbind(name= COG_profile$name,COG_profile[,5:25])
library(janitor)
COG_profile_perc = ns_to_percents(COG_profile_perc)
COG_profile_perc_long <- melt(COG_profile_perc, id=c("name")) 


PROF2 = ggplot(COG_profile_perc_long, aes(x=name,y=value)) + 
    geom_col(aes(fill=variable)) + 
    xlab("Genome") +
    ylab("COG profile") + scale_fill_manual(values = getPalette(21)) +
    theme_classic() + scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()
PROF2


SLCTANV2 = ANVIO[grep('methyl', ANVIO$COG_FUNCTION), ]
head(ANVIO[grep('MTAISETPILDVPILDVKGLAKQFSLHEQNEIVPSCSGVNMLAFPGQLAALTGPTGAGKS', ANVIO$aa_sequence), ])
ANVIO[grep('COG3221', ANVIO$COG_FUNCTION_ACC), ]


SLCTANV2 = ANVIO[grep('ABC', ANVIO$COG_FUNCTION), ]
subselection = data.frame(cbind(as.character(SLCTANV$COG_FUNCTION),as.character(SLCTANV$genome_name)))
colnames(subselection)=c('COG_FUNCTION','name')
dunno = data.frame((table(subselection)))
COG_HEAT = reshape(dunno, idvar = "name", timevar = "COG_FUNCTION", direction = "wide")
rownames(COG_HEAT) =  COG_HEAT$name
COG_HEAT = COG_HEAT[,2:ncol(COG_HEAT)]
COG_HEAT2 = COG_HEAT
COG_HEAT2[COG_HEAT2>1] <-1
COG_HEAT2_ordered <- COG_HEAT2[new_order,]
heatmap.2(as.matrix(t(COG_HEAT2_ordered)),trace="none",col =c('white',"black"),Colv= rep_tree_d,margins = c(1000/ncol(COG_HEAT2), 10),dendrogram="col",cexCol=0.5,cexRow=10/ncol(COG_HEAT2)*3,key=FALSE)



#######################
library(ape)
library(adephylo)
library(phylobase)

ANVIOtree<-read.tree("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/ANVIO_RAxML/RAxML_bipartitions.autosubst100")
ANVIODistMatrix<-cophenetic(ANVIOtree)

ANVIODistMatrix  = ANVIODistMatrix[ order(row.names(ANVIODistMatrix)), ]
ANVIODistMatrix  = ANVIODistMatrix[ , order(colnames(ANVIODistMatrix))]

#distTips(tree)

m2 <- melt(as.matrix(shelldist))[melt(upper.tri(as.matrix(shelldist)))$value,]
m3 <- melt(as.matrix(ANVIODistMatrix))[melt(upper.tri(as.matrix(ANVIODistMatrix)))$value,]

constructedMatrix = data.frame(cbind(m2$value,m3$value))
colnames(constructedMatrix) = c("functional","tree_based")

ggplot(constructedMatrix, aes(functional, tree_based)) + geom_point() + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) 

mantel.rtest(as.dist(ANVIODistMatrix), shelldist, nrepet = 999)

plot(ANVIOnet, "2D")
ANVIOnet = neighborNet(ANVIODistMatrix)

#######################


#---- PLOT TREES ------------------------------
 
my.tree <- panTree(AnvioData)
Anviotree <- panTree(AnvioData,nboot=100)
Pfamtree <- panTree(PfamData,nboot=100)

par(mfrow=c(1,3)) 
plot(tree2)
plot(Anviotree,xlab="PCs Mantattan distance")
plot(Pfamtree,xlab="Pfam domains Mantattan distance")


w <- geneWeights(AnvioData,type="shell")
v <- geneWeights(PfamData,type="shell")

wAnviotree = panTree(AnvioData, scale=0.1, weights=w,nboot = 100)
wPfamtree = panTree(PfamData, scale=0.1, weights=v,nboot = 100)

par(mfrow=c(1,3)) 
plot(tree2)
plot(wAnviotree,xlab="weighted PCs Mantattan distance")
plot(wPfamtree,xlab="weighted Pfam domains Mantattan distance")



#---- BINOMIAL MIXTURE MODEL ------------------------------
# Fitting a binomial mixture model will also produce a conservative estimate of pan-genome size, as well as an estimate of the core size
binomix <- binomixEstimate(AnvioData, K.range=2:11)
binomix2 <- binomixEstimate(PfamData, K.range=2:11)

# Displaying the BIC.table

binomix$BIC.table
binomix2$BIC.table
# Summary of model

summary(binomix)
summary(binomix2)

# Plotting binomial mixture

plot(binomix)

##########################################################################
#						Evolutionary distance ANIb, ANIm, tetra
##########################################################################
ANIb = as.matrix(read.table('~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/JSpecies_ANIb.txt',header=TRUE))
ANIm = as.matrix(read.table('~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/JSpecies_ANIm.txt',header=TRUE))
Tetra = read.table('~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/JSpecies_Tetra.txt',header=TRUE)

CompareM =read.table("~/DATA/MarinobacterGenomics/miscl/CompareM/aai_summary.tsv",header=TRUE)

dat = CompareM[,c("Genome_A","Genome_B","Mean_AAI")]
g <- graph.data.frame(dat, directed=FALSE)
AAI = get.adjacency(g, attr="Mean_AAI", sparse=FALSE)
AAI[AAI<1] <- 100
AAI = AAI[grep("Marinobacter_", rownames(AAI)),grep("Marinobacter_", colnames(AAI))]



DistanceAAI = AAI
DistanceANIb =(ANIb)
DistanceANIm =(ANIm)
DistanceTetra =(Tetra)
PfamDist=dist(PfamData)
PanGenomeDist = dist(AnvioData)

DistanceAAI = DistanceAAI[ order(row.names(DistanceAAI)), ]
DistanceANIb  = DistanceANIb[ order(row.names(DistanceANIb)), ]
DistanceANIm  = DistanceANIm[ order(row.names(DistanceANIm)), ]
DistanceTetra  = DistanceTetra[ order(row.names(DistanceTetra)), ]

DistanceAAI  = DistanceAAI[ , order(colnames(DistanceAAI))]
DistanceANIb  = DistanceANIb[ , order(colnames(DistanceANIb))]
DistanceANIm  = DistanceANIm[ , order(colnames(DistanceANIm))]
DistanceTetra  = DistanceTetra[ , order(colnames(DistanceTetra))]
ANVIODistMatrix  = ANVIODistMatrix[ order(row.names(ANVIODistMatrix)), ]
ANVIODistMatrix  = ANVIODistMatrix[ , order(colnames(ANVIODistMatrix))]
#shelldist  = shelldist[ order(row.names(shelldist)), ]
#shelldist  = shelldist[ , order(colnames(shelldist))]
#clouddist  = shelldist[ order(row.names(clouddist)), ]
#clouddist  = shelldist[ , order(colnames(clouddist))]
PfamDist  = PfamDist[ order(row.names(PfamDist)), ]
PfamDist  = PfamDist[ , order(colnames(PfamDist))]
PanGenomeDist  = PanGenomeDist[ order(row.names(PanGenomeDist)), ]
PanGenomeDist  = PanGenomeDist[ , order(colnames(PanGenomeDist))]

DistanceAAI = as.matrix(DistanceAAI)
DistanceANIb = as.matrix(DistanceANIb)
DistanceANIm = as.matrix(DistanceANIm)
DistanceTetra = as.matrix(DistanceTetra)

m2 <- melt(DistanceANIb)[melt(upper.tri(DistanceANIb))$value,]
m3 <- melt(DistanceTetra)[melt(upper.tri(DistanceTetra))$value,]
m4 <- melt(DistanceANIm)[melt(upper.tri(DistanceANIm))$value,]
m5 <- melt(as.matrix(ANVIODistMatrix))[melt(upper.tri(as.matrix(ANVIODistMatrix)))$value,]
m6 <- melt(as.matrix(shelldist))[melt(upper.tri(as.matrix(shelldist)))$value,]
m7 <- melt(as.matrix(clouddist))[melt(upper.tri(as.matrix(clouddist)))$value,]
m8 <- melt(as.matrix(PanGenomeDist))[melt(upper.tri(as.matrix(PanGenomeDist)))$value,]
m9 <- melt(as.matrix(PfamDist))[melt(upper.tri(as.matrix(PfamDist)))$value,]
m10 <- melt(as.matrix(DistanceAAI))[melt(upper.tri(as.matrix(DistanceAAI)))$value,]

m11 <- melt(as.matrix(GeoDist))[melt(upper.tri(as.matrix(GeoDist)))$value,]
m2 <- melt(DistanceANIb)[melt(upper.tri(DistanceANIb))$value,]

testco = data.frame(cbind('ANIb' = m2$value,"Geo"=m11$value))

constructedMatrix = data.frame(cbind('ANIb' = m2$value, 'Tetra'=m3$value,'ANIm'=m4$value,'ANVIO'=m5$value,'shell'=m6$value,'cloud'=m7$value,'PanGenome'=m8$value,'Pfam'=m9$value,"AAI"=m10$value))

test = melt(constructedMatrix,id="ANIb")

a0 = ggplot(constructedMatrix, aes(ANIb, AAI)) + geom_point() + geom_line(data = log.model.df, aes(x, y, color = "Log Model"), size = 1, linetype = 2) + scale_x_reverse()+ scale_y_reverse()+theme_classic()
a1 = ggplot(constructedMatrix, aes(ANIb, Tetra)) + geom_point() + geom_smooth(method="lm", formula= (y ~ exp(x)), se=FALSE,color="darkgreen") + scale_x_reverse()+ scale_y_reverse()+theme_classic()
a2 = ggplot(constructedMatrix, aes(ANIb, ANVIO)) + geom_point()  + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="darkgreen", se=F, fullrange=T) + scale_x_reverse()+ theme_classic()
a3 = ggplot(constructedMatrix, aes(ANIb, PanGenome)) + geom_point() + geom_smooth(method="lm", color="darkgreen", se=F, fullrange=T) +scale_x_reverse()+  scale_y_reverse()+theme_classic()
a4 = ggplot(constructedMatrix, aes(ANIb, Pfam)) + geom_point() +  geom_smooth(method="lm", color="darkgreen", se=F, fullrange=T) +scale_x_reverse()+ scale_y_reverse()+ theme_classic()

multiplot(a0,a1,a2,a3,a4,cols=5)



m2 <- melt(DistanceANIb)[melt(upper.tri(DistanceANIb))$value,]
m3 <- melt(DistanceTetra)[melt(upper.tri(DistanceTetra))$value,]
m4 <- melt(scale(DistanceANIm))[melt(upper.tri(scale(DistanceANIm)))$value,]
m5 <- melt(as.matrix(scale(ANVIODistMatrix)))[melt(upper.tri(as.matrix(scale(ANVIODistMatrix))))$value,]
m6 <- melt(as.matrix(scale(shelldist)))[melt(upper.tri(as.matrix(scale(shelldist))))$value,]
m7 <- melt(as.matrix(scale(clouddist)))[melt(upper.tri(as.matrix(scale(clouddist))))$value,]

ggplot(test, aes(ANIb, value,color=variable)) + geom_point() + facet_wrap(~variable,scale='free')+scale_x_reverse()

test2 = melt(constructedMatrix,id="ANVIO")

ggplot(test2, aes(ANVIO, value,color=variable)) + geom_point() + facet_wrap(~variable,scale='free')+scale_x_reverse()


layout3d(matrix(1:4, 2, 2),sharedMouse = TRUE)
    
plot3d(constructedMatrix[,1:3], size=15, type='p')
grid3d(c("x", "y+", "z"))
next3d()

plot3d(constructedMatrix[,2:4], col= 'black', size=15, type='p')
grid3d(c("x", "y+", "z"))
next3d()

plot3d(constructedMatrix[,3:5], col= 'black', size=15, type='p')
grid3d(c("x", "y+", "z"))
next3d()

plot3d(constructedMatrix[,4:6], col= 'black', size=15, type='p')
grid3d(c("x", "y+", "z"))
next3d()

plot3d(cbind(constructedMatrix[,1],constructedMatrix[,5:6]), col= 'black', size=15, type='p')
grid3d(c("x", "y+", "z"))
next3d()

multiplot(a1,a2,a3,cols=3)

par(mfrow=c(2,2))
h1 = heatmap.2(ANIb,trace="none",scale="none",col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))
h2 = heatmap.2(ANIm,trace="none",scale="none",col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))
h4 = heatmap.2(as.matrix(AAI),trace="none",scale="none",col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))
h3 = heatmap.2(as.matrix(Tetra),trace="none",scale="none",col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))

multiplot(h1,h2,h3,h4,cols=4)

h1 = heatmap.2(ANIb,trace="none",scale="none",hclustfun = function(x) hclust(x,method = 'ward.D2'),col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))
h2 = heatmap.2(ANIm,trace="none",scale="none",hclustfun = function(x) hclust(x,method = 'ward.D2'),col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))
h3 = heatmap.2(as.matrix(Tetra),trace="none",scale="none",hclustfun = function(x) hclust(x,method = 'ward.D2'),col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))
h4 = heatmap.2(as.matrix(AAI),trace="none",scale="none",hclustfun = function(x) hclust(x,method = 'ward.D2'),col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))

multiplot(h1,h2,h3,cols=3)

h1 = heatmap.2(ANIb,trace="none",scale="none",hclustfun = function(x) hclust(x,method = 'ward.D'),col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))
h2 = heatmap.2(ANIm,trace="none",scale="none",hclustfun = function(x) hclust(x,method = 'ward.D'),col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))
h3 = heatmap.2(as.matrix(Tetra),trace="none",scale="none",hclustfun = function(x) hclust(x,method = 'ward.D'),col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))
multiplot(h1,h2,h3,cols=3)

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

# perform clustering on rows and columns
cl.row <- hclustfunc(distfunc(ANIb))
cl.col <- hclustfunc(distfunc(t(ANIb)))

# extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
gr.row <- cutree(cl.row, 9)
gr.col <- cutree(cl.col, 9)

# require(RColorBrewer)
col1 <- brewer.pal(8, "Set1")
col2 <- brewer.pal(8, "Pastel1")

h1 = heatmap.2(ANIb,trace="none",scale="none",hclustfun = function(x) hclust(x,method = 'ward.D'),col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15),RowSideColors=pcacolors, ColSideColors=pcacolors)



mantel.rtest(as.dist(ANIb), as.dist(ANIm), nrepet = 999)
mantel.rtest(as.dist(ANIb),as.dist(Tetra), nrepet = 999)
mantel.rtest(as.dist(ANIm), as.dist(Tetra), nrepet = 999)
mantel.rtest(as.dist(ANIb), as.dist(AAI), nrepet = 999)


#convert to long format
longANIb = melt(ANIb)

#ANI treshold of 96 (95) is considered to deliniate species 
Treshold_longANIb = longANIb[(longANIb $value>95.5),]
Treshold_longANIb = Treshold_longANIb[(Treshold_longANIb$value<100),]

Treshold_longANIb = longANIb[(longANIb $value<71),]
Treshold_longANIb = Treshold_longANIb[(Treshold_longANIb$value<100),]

STW2_longANIb = longANIb[(longANIb$Var1=='Marinobacter_hydrocarbonoclasticus_STW2'),]



ggplot(constructedMatrix, aes(functional, tree_based)) + geom_point() + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) 

mantel.rtest(as.dist(ANVIODistMatrix), shelldist, nrepet = 999)







##########################################################################
#						COMPARE ANVIO AND CHECKM
##########################################################################
library(RColorBrewer)

CheckM = read.table("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/QC/qa2Marinobacter.txt",header=TRUE)
ANVIOcheck = read.table("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/QC/MARINOBACTER-samples-information.txt",header=TRUE)
BUSCOcheck = read.table("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/QC/Marinobacter_BUSCO_SUMMARY.txt",header = TRUE)
CompareCharacter = read.table("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/QC/CheckM_vs_others.txt",header = TRUE,sep="\t")



CheckM2 = subset(CheckM, Completeness > 98)
CheckM2 = CheckM2[order(CheckM2$Bin_Id),]

BUSCOcheck = BUSCOcheck[order(BUSCOcheck$genome),]
BUSCOcheck = cbind(BUSCOcheck, Bcompleteness = (148-BUSCOcheck$M)/148*100)
BUSCOcheck = cbind(BUSCOcheck,C=BUSCOcheck$S + BUSCOcheck$D)

ANVIOcheck = ANVIOcheck[order(ANVIOcheck$samples),]

multiData = cbind(CheckM$Genome_size,CheckM$GC,CheckM$predicted_genes)
cols <- palette(brewer.pal(11,'Set1'))[CheckM$Clade]
rownames(multiData) = CheckM$Bin_Id
colnames(multiData) = c("Genome_size","GC","Predicted_genes")

CheckM2 = CheckM2[CheckM2$Bin_Id != "Marinobacter_dagiaonensis_CGMCC_1_9167",]
BUSCOcheck = BUSCOcheck[BUSCOcheck$genome != "Marinobacter_dagiaonensis_CGMCC_1.9167",]

plot(CheckM2$Completeness,ANVIOcheck$percent_complete)
plot(CheckM2$Genome_size,ANVIOcheck$total_length)
plot(CheckM2$GC,ANVIOcheck$gc_content)
plot(CheckM2$predicted_genes,ANVIOcheck$num_genes)
plot(CheckM2$coding_density,ANVIOcheck$num_genes_per_kb)
combinedCheck = cbind(CheckM2,ANVIOcheck,BUSCOcheck)

p1 <- ggplot(combinedCheck, aes(x= Genome_size/1000000, y=total_length/1000000)) + xlab("Genome size (Mbp) CheckM")+ ylab("Genome size (Mbp) Anvio") + geom_point() + geom_smooth(method="lm",color = "orange", size = 0.75)+theme_classic()

p2 <- ggplot(combinedCheck, aes(x= GC, y=gc_content*100)) + geom_point() + xlab("%GC CheckM")+ ylab("%GC Anvio") + geom_smooth(method="lm",color = "orange", size = 0.75)+theme_classic()
p3 <- ggplot(combinedCheck, aes(x= predicted_genes, y= num_genes)) + geom_point() + xlab("predicted genes - CheckM")+ ylab("predicted genes - Anvio") + geom_smooth(method="lm",color = "orange", size = 0.75)+theme_classic()

p4 <- ggplot(combinedCheck, aes(x= coding_density/100, y= num_genes_per_kb)) + geom_point() + xlab("coding density - CheckM")+ ylab("coding density - Anvio") + geom_smooth(method="lm",color = "orange", size = 0.75)+theme_classic()

p5 <- ggplot(combinedCheck, aes(x= Completeness, y=percent_complete)) + xlab("%Completeness - CheckM")+ ylab("%Completeness - Anvio") + geom_point() + geom_smooth(method="lm",color = "orange", size = 0.75)+theme_classic()

p6 <- ggplot(combinedCheck, aes(x= Completeness, y=Bcompleteness)) + xlab("%Completeness - CheckM")+ ylab("%Completeness - BUSCO") + geom_point() + geom_smooth(method="lm",color = "orange", size = 0.75)+theme_classic()


rcorr(combinedCheck $Genome_size, combinedCheck $total_length, type="pearson")
summary(lm(combinedCheck $total_length ~ combinedCheck $Genome_size))
summary(lm(combinedCheck $gc_content*100 ~ combinedCheck $GC))
summary(lm(combinedCheck $num_genes ~ combinedCheck $predicted_genes))
summary(lm(combinedCheck$num_genes_per_kb ~ combinedCheck$coding_density))
summary(lm(combinedCheck $percent_complete ~ combinedCheck $Completeness))
summary(lm(combinedCheck $Bcompleteness ~ combinedCheck $Completeness))

multiplot(p1, p4, p2, p5,p3,p6, cols=3)
multiplot(p1, p2, p3, cols=3)
multiplot(p4, p5, cols=3)


############################################################
# Compare the CheckM obtained statistics with the online data
###############################################################


p1 <- ggplot(CompareCharacter, aes(x= Genome_size/1000000, y= Genome_Length/1000000)) + xlab("Genome size (Mbp) CheckM")+ ylab("Genome size (Mbp) NCBI") + geom_point() + geom_smooth(method="lm",color = "orange", size = 0.75)+theme_classic()

p2 <- ggplot(CompareCharacter, aes(x= GC, y=GC_Content*100)) + geom_point() + xlab("%GC CheckM")+ ylab("%GC NCBI") + geom_smooth(method="lm",color = "orange", size = 0.75)+theme_classic()
p3 <- ggplot(CompareCharacter, aes(x= predicted_genes, y= PATRIC_CDS)) + geom_point() + xlab("predicted genes - CheckM")+ ylab("predicted genes - PATRIC") + geom_smooth(method="lm",color = "orange", size = 0.75)+theme_classic()
p4 <- ggplot(CompareCharacter, aes(x= predicted_genes, y= RefSeq_CDS)) + geom_point() + xlab("predicted genes - CheckM")+ ylab("predicted genes - RefSeq") + geom_smooth(method="lm",color = "orange", size = 0.75)+theme_classic()
multiplot(p1, p2, p3,p4, cols=2)

rcorr(CompareCharacter $Genome_size, CompareCharacter $Genome_Length, type="pearson")

summary(lm(CompareCharacter $Genome_size ~ CompareCharacter $Genome_Length))
summary(lm(CompareCharacter $GC*100 ~ CompareCharacter $GC_Content))
summary(lm(CompareCharacter $predicted_genes ~ CompareCharacter $PATRIC_CDS))
summary(lm(CompareCharacter $predicted_genes ~ CompareCharacter $RefSeq_CDS))


######## import the phylogenetic tree to sort the CheckM values

RAxMLANVIO = read.tree("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/ANVIO_RAxML/RAxML_bipartitions.autosubst100")
RAxMLANVIORooted =  midpoint.root(RAxMLANVIO)

tree2 <- ladderize(RAxMLANVIORooted, right = FALSE)


######## MIDPOINT ROOT THE TREE
RAxMLANVIORooted =  midpoint.root(RAxMLANVIO)

######## sort branches based on lenght, increasing
tree2 <- ladderize(RAxMLANVIORooted, right = FALSE)
tr = plot(tree2)

######## Obtain the ordering of the tip labels
is_tip <- tree2$edge[,2] <= length(tree2$tip.label)
ordered_tips <- tree2$edge[is_tip, 2]
tree2$tip.label[ordered_tips]

######## plotting the different CheckMvalues using the tip label order
q = ggplot(CheckM2, aes(x=Bin_Id, y=GC)) + 
    geom_point() +      # Thinner lines
    xlab("Genome") +
    ylab("GC content") +
    theme_bw() + scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()

q1 = ggplot(CheckM, aes(x=Bin_Id, y=GC)) + 
    geom_point() +      # Thinner lines
    xlab("Genome") +
    ylab("GC content") +theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()

q2  = ggplot(CheckM, aes(x=Bin_Id, y=0.000001*Genome_size)) + 
    geom_point() +      # Thinner lines
    ylab("Genome size (Mbp)") + theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()

q3  = ggplot(CheckM, aes(x=Bin_Id, y=predicted_genes)) + 
    geom_point() +      # Thinner lines
    ylab("# predicted genes") +  theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()
    
q4  = ggplot(CheckM, aes(x=Bin_Id, y= coding_density)) + 
    geom_point() +      # Thinner lines
    ylab("coding density") + theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()
    

q7  = ggplot(CheckM, aes(x=Bin_Id, y= contigs)) + 
    geom_point() +      # Thinner lines
    ylab("# contigs") + theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()

q8  = ggplot(CheckM, aes(x=Bin_Id, y= Completeness)) + 
    geom_point() +      # Thinner lines
    ylab("Completeness (%)") + theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + ylim(98,100) + 
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()

p1 = ggplot(DF1, aes(x = name, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + xlab("Genome") + ylab("number of unique PCs") + coord_flip() + scale_x_discrete(limits=tree2$tip.label[ordered_tips]) + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +theme(legend.position="none")+theme(axis.line = element_line(color="black", size = 0.5))

DF2$name = DF1$name
p2 = ggplot(DF2, aes(x = name, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + xlab("Genome") + ylab("number of unique PCs") + coord_flip() + scale_x_discrete(limits=tree2$tip.label[ordered_tips])+scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) +  theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + theme(legend.position="none") +theme(axis.line = element_line(color="black", size = 0.5))




multiplot(tr,p1,p2,q1,q2,q3,q4,q8,cols=8)


multiplot(q1,q2,q3,q4,q8,cols=5)




q = ggplot(CheckM2, aes(x=Bin_Id, y=GC)) + 
    geom_col() +      # Thinner lines
    xlab("Genome") +
    ylab("GC content") +
    theme_bw() + scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()

q1 = ggplot(CheckM, aes(x=Bin_Id, y=GC)) + 
    geom_col() +      # Thinner lines
    xlab("Genome") +
    ylab("GC content") +theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()

q2  = ggplot(CheckM, aes(x=Bin_Id, y=0.000001*Genome_size)) + 
    geom_col() +      # Thinner lines
    ylab("Genome size (Mbp)") + theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()

q3  = ggplot(CheckM, aes(x=Bin_Id, y=predicted_genes)) + 
    geom_col() +      # Thinner lines
    ylab("# predicted genes") +  theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()
    
q4  = ggplot(CheckM, aes(x=Bin_Id, y= coding_density)) + 
    geom_col() +      # Thinner lines
    ylab("coding density") + theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()
    

q7  = ggplot(CheckM, aes(x=Bin_Id, y= contigs)) + 
    geom_col() +      # Thinner lines
    ylab("# contigs") + theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()

q8  = ggplot(CheckM, aes(x=Bin_Id, y= Completeness)) + 
    geom_col() +      # Thinner lines
    ylab("Completeness (%)") + theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + ylim(98,100) + 
     scale_x_discrete(limits=tree2$tip.label[ordered_tips]) +  coord_flip()

p1 = ggplot(DF1, aes(x = name, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + xlab("Genome") + ylab("number of unique PCs") + coord_flip() + scale_x_discrete(limits=tree2$tip.label[ordered_tips]) + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +theme(legend.position="none")+theme(axis.line = element_line(color="black", size = 0.5))

DF2$name = DF1$name
p2 = ggplot(DF2, aes(x = name, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + xlab("Genome") + ylab("number of unique PCs") + coord_flip() + scale_x_discrete(limits=tree2$tip.label[ordered_tips])+scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) +  theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + theme(legend.position="none") +theme(axis.line = element_line(color="black", size = 0.5))



#######################
#		STAMPomaniac
#######################
library(dplyr)

STAMP1 = read.table(file = '~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PC_algicola-clade_FDR.tsv', sep = '\t', header = TRUE)
STAMP2 = read.table(file = '~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PC_psycrhophylus-clade_FDR.tsv', sep = '\t', header = TRUE)
STAMP1 = read.table(file = '~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PC_adhaerens-clade_FDR.tsv', sep = '\t', header = TRUE)
STAMP1 = read.table(file = '~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PC_antarcticus-clade_FDR.tsv', sep = '\t', header = TRUE)
STAMP1 = read.table(file = '~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PC_excellens-clade_FDR.tsv', sep = '\t', header = TRUE)
STAMP1 = read.table(file = '~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PC_hydrocarbonoclasticus-clade_FDR.tsv', sep = '\t', header = TRUE)
STAMP1 = read.table(file = '~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PC_lipolyticus-clade_FDR.tsv', sep = '\t', header = TRUE)
STAMP1 = read.table(file = '~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PC_viniformus-clade_FDR.tsv', sep = '\t', header = TRUE)

ANVIO <- read.table('~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/MARINOBACTER_protein_clusters_summary.txt',sep='\t',header=TRUE)

DEDUPLICATED =  subset(ANVIO, !duplicated(protein_cluster_id))
colnames(DEDUPLICATED)=c("unique_id","gr","bin_name","genome_name","gene_callers_id","COG_CATEGORY_ACC","COG_CATEGORY","COG_FUNCTION_ACC","COG_FUNCTION","aa_sequence" )

JOINED <- right_join(DEDUPLICATED, STAMP1, by=c("gr"))
JOINED = JOINED[order(-JOINED$X95.0..upper.CI),]


STAMP2 = read.table(file = '~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PC_algicola-clade_FDR.tsv', sep = '\t', header = TRUE)

STAMP2 = read.table(file = '~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/PC_psychrophylus-clade_FDR.tsv', sep = '\t', header = TRUE)


JOINED <- right_join(DEDUPLICATED, STAMP2, by=c("gr"))
JOINED = JOINED[order(-JOINED$X95.0..upper.CI),]


JOINED[grep('C', JOINED$COG_CATEGORY_ACC),]
write.csv(JOINED[grep('C', JOINED$COG_CATEGORY_ACC),],'~/Algicola_C_enriched_list.csv')


nrCOGS = c()
COG_list = c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W" )
for(i in COG_list){
	#nrow(JOINED[grep(i, JOINED$COG_CATEGORY_ACC),])
	nrCOGS = c(nrCOGS,nrow(JOINED[grep(i, JOINED$COG_CATEGORY_ACC),]))
	write.csv(JOINED[grep(i, JOINED$COG_CATEGORY_ACC),c(2,6:19, seq(21, 200, by = 3))],paste('~/Psychro_' , i , '_enriched_list.csv',sep=''))

}
names(nrCOGS)=COG_list

COG_table2 = data.frame(cbind(COG_CATEGORY=as.character((JOINED$COG_CATEGORY)),group = rep('vinifirmus',nrow(JOINED))))
COG_ultimate= rbind(COG_ultimate, COG_table2)

ggplot(COG_ultimate) + geom_bar(aes(x=COG_CATEGORY,fill=group),stat = "count") + theme(axis.text.x=element_text(angle = 90, hjust = 0))

 ggplot(COG_ultimate) + geom_bar(aes(x=COG_CATEGORY,fill=group),stat = "count",position="dodge") + theme(axis.text.x=element_text(angle = 90, hjust = 0))

SELECT = cbind(JOINED[,2],JOINED[,6:10],group = JOINED[,11], others = JOINED[,13])

SELECT_molten = melt(head(SELECT,75),ids=c("unique_id","gr","bin_name","genome_name","gene_callers_id","COG_CATEGORY_ACC","COG_CATEGORY","COG_FUNCTION_ACC","COG_FUNCTION","aa_sequence","algicola..std..dev.....","X.All.other.samples...std..dev....." ))

p2 = ggplot(SELECT_molten, aes(x = COG_FUNCTION, y = value,fill=variable)) + 
  geom_bar(stat = "identity",position="dodge") +  theme(axis.text.y  = element_text(size=7)) + coord_flip()
p2  
  
SELECT_molten = melt(SELECT,ids=c("unique_id","gr","bin_name","genome_name","gene_callers_id","COG_CATEGORY_ACC","COG_CATEGORY","COG_FUNCTION_ACC","COG_FUNCTION","aa_sequence","algicola..std..dev.....","X.All.other.samples...std..dev....." ))

p2 = ggplot(SELECT_molten, aes(x = COG_FUNCTION, y = value,fill=variable)) + 
  geom_bar(stat = "identity",position="dodge") +  theme(axis.text.y  = element_text(size=1)) + coord_flip()
p2  
  
  

c("unique_id","gr","bin_name","genome_name","gene_callers_id","COG_CATEGORY_ACC","COG_CATEGORY","COG_FUNCTION_ACC","COG_FUNCTION","aa_sequence","algicola..std..dev.....","X.All.other.samples...std..dev....." )

as.character(unique(ANVIO[ANVIO$protein_cluster_id %in% STAMP1[,1],9]))

ann = data.frame(cbind(as.character(head(STAMP1[,1])),annotation=as.character(unique(ANVIO[ANVIO$protein_cluster_id %in% as.character(head(STAMP1[,1])),9])),head(STAMP1[,2]),head(STAMP1[,3])))

p2 = ggplot(ann, aes(x = annotation, y = V3)) + geom_point() + coord_flip()


 p2 = ggplot(ann, aes(x = ann, y = V3)) + 
  geom_bar(stat = "identity",position="dodge") + coord_flip()


cbind(STAMP1[,1],STAMP1[,2],STAMP1[,4])
plot(STAMP1[,2],STAMP1[,4])

Stamp_molten = melt(data.frame(cbind(as.character(STAMP1[,1]),STAMP1[,2],STAMP1[,4])),id=c("X1"))
Stamp_molten = melt(data.frame(cbind(as.character(head(STAMP1[,1],100)),head(STAMP1[,2],100),head(STAMP1[,4],100))),id=c("X1"))

Stamp_molten = melt(data.frame(cbind(as.character(tail(STAMP1[,1],50)),tail(STAMP1[,2],50),tail(STAMP1[,4],50))),id=c("X1"))

cbind(Stamp_molten,as.character(unique(ANVIO[ANVIO$protein_cluster_id %in% Stamp_molten$X1,1]))


as.character(unique(ANVIO[ANVIO$protein_cluster_id=="PC_00005024",9]))




Stamp_molten


 p2 = ggplot(Stamp_molten, aes(x = X1, y = as.numeric(paste(value)), fill = variable)) + 
  geom_bar(stat = "identity",position="dodge") + coord_flip()

p2 = ggplot(Stamp_molten, aes(x = name, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + xlab("Genome") + ylab("number of unique PCs") + coord_flip() + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) +  theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + theme(legend.position="none") +theme(axis.line = element_line(color="black", size = 0.5))


#######################
#		EGSEA
#######################
source("http://bioconductor.org/biocLite.R")
biocLite('EGSEA')
biocLite(c("PADOG", "GSVA", "AnnotationDbi", "topGO", "pathview",
    "gage", "globaltest", "limma", "edgeR", "safe", "org.Hs.eg.db",
    "org.Mm.eg.db", "org.Rn.eg.db"))
biocLite("EGSEAdata")
install_bitbucket("malhamdoosh/egseadata", ref = "Stable_Release")
install_bitbucket("malhamdoosh/egsea", ref = "Devel_Release")

library(EGSEA) library(EGSEAdata)
data(il13.data)
v = il13.data$voom
names(v)






#######################
#		rBLAST
#######################
source("https://bioconductor.org/biocLite.R")
biocLite()
library(devtools)
install_github("mhahsler/rBLAST")
Library(rBLAST)
library(seqinr)

download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz","16SMicrobial.tar.gz", mode='wb')
untar("16SMicrobial.tar.gz", exdir="16SMicrobialDB")
seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",package="rBLAST"))
bl <- blast(db="./16SMicrobialDB/16SMicrobial")
cl <- predict(bl, seq[1,])


seq_1_DNAbin <- read.GenBank("NR_025644")
attr(seq_1_DNAbin, "species")
str(seq_1_DNAbin)
lizards_accession_numbers <- c("JF806202", "HM161150", "FJ356743", "JF806205", 
                               "JQ073190", "GU457971", "FJ356741", "JF806207",
                               "JF806210", "AY662592", "AY662591", "FJ356748",       
                               "JN112660", "AY662594", "JN112661", "HQ876437", 
                               "HQ876434", "AY662590", "FJ356740", "JF806214", 
                               "JQ073188", "FJ356749", "JQ073189", "JF806216",
                               "AY662598", "JN112653", "JF806204", "FJ356747",
                               "FJ356744", "HQ876440", "JN112651", "JF806215",
                               "JF806209") 
lizards_sequences <- read.GenBank(lizards_accession_numbers) 

str(lizards_sequences) 
attr(lizards_sequences, "species")

write.dna(lizards_sequences, file ="lizard_fasta_1.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)
lizard_seq_seqinr_format <- read.fasta(file = "lizard_fasta_1.fasta", seqtype = "DNA",as.string = TRUE, forceDNAtolower = FALSE)
write.fasta(sequences = lizard_seq_seqinr_format, names = lizards_sequences_GenBank_IDs,nbchar = 10, file.out = "lizard_seq_seqinr_format.fasta")

install.packages ("rentrez")
library (rentrez)
lizard <- "Basiliscus basiliscus[Organism]" 
#We want a character vector!
lizard_search <- entrez_search(db="nuccore", term=lizard, retmax=40)
lizard_search
lizard_search$ids #gives you the NCBI ids
#gets your sequences as a character vector
lizard_seqs <- entrez_fetch(db="nuccore", id=lizard_search$ids, rettype="fasta")
lizard_seqs

lizards_char = as.character(lizards_sequences)


marinobact <- "(Marinobacter[Organism]) AND 16S" 
#We want a character vector!
marinobact_search <- entrez_search(db="nuccore", term= marinobact, retmax=40)
marinobact_search
marinobact_search $ids #gives you the NCBI ids
#gets your sequences as a character vector
marinobact_seqs <- entrez_fetch(db="nuccore", id= marinobact_search$ids, rettype="fasta")
marinobact_seqs

lizards_sequences <- read.GenBank(lizards_accession_numbers) 



data(Laurasiatherian)
class(Laurasiatherian)
LauraChar <- as.character(Laurasiatherian)
Laura <- phyDat(LauraChar, return.index=TRUE)

mammals = read.dna("~/lizard_fasta_1.fasta", format="fasta")
mammals_phyDat <- phyDat(mammals, type = "DNA", levels = NULL)
mammals10 <- subset(mammals_phyDat, 1:10)
mammals10_phyDat <- phyDat(mammals10, type = "DNA", levels = NULL)


modelTest(Laura)
dna_dist <- dist.ml(Laura, model="JC")
upgma = upgma(dna_dist)
nj = nj(Laura)
fit <- pml(nj, Laura)

fitJC <- optim.pml(fit, model = "JC", rearrangement ="stochastic")
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")


write.tree(bs, file="bootstrap_example.tre")

pheatmap((table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% c('COG1830', 'COG0837', 'COG0205', 'COG0191', 'COG0166', 'COG0057', 'COG0149','COG0126' , 'COG0148','COG0469','COG0588','COG0696','COG3635'))[,c('COG_FUNCTION','genome_name')]))),cex=0.6)

pheatmap((table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% c('COG0538', 'COG0567', 'COG1048', 'COG0473', 'COG0372', 'COG0074', 'COG0045','COG0479' , 'COG1053','COG0508','COG1249','COG0039','COG0114','COG1951','COG1838'))[,c('COG_FUNCTION','genome_name')]))),cex=0.6)

pheatmap((table(droplevels(ANVIO[grep('pyruvate',ANVIO$COG_FUNCTION),c('COG_FUNCTION','genome_name')]))),cex=0.6)

mat = table(droplevels(ANVIO[grep('pyruvate',ANVIO$COG_FUNCTION),c('COG_FUNCTION','genome_name')]))
mat = table(droplevels(ANVIO[grep('transport',ANVIO$COG_FUNCTION),c('COG_FUNCTION_ACC','genome_name')]))

mat = table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% c('COG0538', 'COG0567', 'COG1048', 'COG0473', 'COG0372', 'COG0074', 'COG0045','COG0479' , 'COG1053','COG0508','COG1249','COG0039','COG0114','COG1951','COG1838'))[,c('COG_FUNCTION','genome_name')]))


clade_order <- order.dendrogram(rep_tree_d)
clade_name <- labels(rep_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(AccesoryDF))

mat = t(table(droplevels(ANVIO[grep('ABC',ANVIO$COG_FUNCTION),c('COG_FUNCTION_ACC','genome_name')])))
mat = t(table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% c('COG1830', 'COG0837', 'COG0205', 'COG0191', 'COG0166', 'COG0057', 'COG0149','COG0126' , 'COG0148','COG0469','COG0588','COG0696','COG3635'))[,c('COG_FUNCTION','genome_name')])))

mat = t(table(droplevels(ANVIO[grep('pyruvate',ANVIO$COG_FUNCTION),c('COG_FUNCTION','genome_name')])))
mat = t(table(droplevels(COREmetadata[grep('Q',ANVIO$COG_CATEGORY),c('protein_cluster_id','genome_name')])))



mat_HEAT = mat[,2:ncol(mat)]
mat_HEAT2 = mat_HEAT
mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]

#heatmap.2(mat_HEAT2_ordered,trace="none",col = inferno(75),Rowv= rep_tree_d,labCol = TRUE,margins = c(30, 20),dendrogram="row")
heatmap.2(mat_HEAT2_ordered,trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)


mat = t(table(droplevels(ANVIO[grep('pyruvate',ANVIO$COG_FUNCTION),c('COG_FUNCTION_ACC','genome_name')])))
mat_HEAT = mat[,2:ncol(mat)]
mat_HEAT2 = mat_HEAT
mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]
heatmap.2(mat_HEAT2_ordered,trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)



mat = t(table(droplevels(ANVIO[grep('R',ANVIO$COG_CATEGORY),c('COG_FUNCTION_ACC','genome_name')])))
mat_HEAT = mat[,2:ncol(mat)]
mat_HEAT2 = mat_HEAT
mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]
heatmap.2(mat_HEAT2_ordered,trace="none",col = inferno(16),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)

mat = t(table(droplevels(ANVIO[grep('N',ANVIO$COG_CATEGORY),c('COG_FUNCTION','genome_name')])))
mat_HEAT = mat[,2:ncol(mat)]
mat_HEAT2 = mat_HEAT
mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]
heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>5],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)
heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>5],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5,dendrogram="row",Colv = FALSE)


mat = t(table(droplevels(COREmetadata[grep('X', COREmetadata $COG_CATEGORY),c('COG_FUNCTION','genome_name')])))
mat_HEAT = mat[,2:ncol(mat)]
mat_HEAT2 = mat_HEAT
mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]
heatmap.2(mat_HEAT2_ordered,trace="none",col = inferno(16),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)


mat = t(table(droplevels(ANVIO[,c('COG_FUNCTION','genome_name')])))
mat_HEAT = mat[,2:ncol(mat)]
mat_HEAT2 = mat_HEAT
mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]
heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>5],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)
heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>120],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5,dendrogram="row",Colv = FALSE)

heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>250],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5,dendrogram="row",Colv = FALSE)

heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>1000],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)
heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>750],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)
heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>600],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)

heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>300],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)
heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>600],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)

heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>300],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)

mat = t(table(droplevels(ANVIO[,c('COG_CATEGORY','genome_name')])))


ANVIO[grep('C',ANVIO$COG_CATEGORY),c('COG_FUNCTION','genome_name')]


mat = t(table(droplevels(ANVIO$COG_CATEGORY[grep('', ANVIO$COG_CATEGORY),c('COG_FUNCTION','genome_name')])))
mat_HEAT = mat[,2:ncol(mat)]
mat_HEAT2 = mat_HEAT
mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]
heatmap.2(mat_HEAT2_ordered,trace="none",col = inferno(16),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)


ANVIO2 = ANVIO
ANVIO2$combined = as.factor(paste(ANVIO$COG_FUNCTION,ANVIO$COG_FUNCTION_ACC,sep='___'))
ANVIO2$combined2 = as.factor(paste('>',ANVIO$genome_name,ANVIO$unique_id,sep='-'))
ANVIO2$combined3 = as.factor(paste('>',ANVIO$genome_name,sep=''))

mat = t(table(droplevels(COREmetadata[grep('X', COREmetadata $COG_CATEGORY),c('COG_FUNCTION','genome_name')])))
mat_HEAT = mat[,1:ncol(mat)]
mat_HEAT2 = mat_HEAT
mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]
heatmap.2(mat_HEAT2_ordered,trace="none",col = inferno(16),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)
heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>50],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)

heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>200],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)

heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>500],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)

heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>300],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20, 20),srtCol=45,cexRow=0.5,cexCol=0.5)



mat = t(table(droplevels(subset(ANVIO, COG_FUNCTION  %in% c('transposase', 'integrase','transpos'))[,c('COG_FUNCTION','genome_name')])))
mat = t(table(droplevels(subset(ANVIO, COG_FUNCTION  %in% c('Transposase', 'Integrase','transpos','integr'))[,c('COG_FUNCTION_ACC','genome_name')])))

mat = t(table(droplevels(subset(ANVIO2, COG_FUNCTION  %in% c('Transposase', 'Integrase','transpos','integr'))[,c('combined','genome_name')])))

mat = t(table(droplevels(ANVIO2[grep('V', ANVIO2$COG_CATEGORY),c('combined','genome_name')])))
mat = t(table(droplevels(ANVIO2[grep('L', ANVIO2$COG_CATEGORY),c('combined','genome_name')])))
mat = t(table(droplevels(ANVIO2[grep('Q', ANVIO2$COG_CATEGORY),c('combined','genome_name')])))
mat = t(table(droplevels(ANVIO2[grep('Q', ANVIO2$COG_CATEGORY),c('combined','genome_name')])))
mat = t(table(droplevels(ANVIO2[grep('E', ANVIO2$COG_CATEGORY),c('combined','genome_name')])))
mat = t(table(droplevels(ANVIO2[grep('K', ANVIO2$COG_CATEGORY),c('combined','genome_name')])))

mat = t(table(droplevels(ANVIO2[,c('combined','genome_name')])))

mat = t(table(droplevels(ANVIO2[grep('V', ANVIO2$COG_CATEGORY),c('combined','genome_name')])))



COG_list = c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W" )
for(i in COG_list){
	mat = t(table(droplevels(ANVIO2[grep(i, ANVIO2$COG_CATEGORY),c('combined','genome_name')])))
	mat_HEAT = mat[,2:ncol(mat)]
	mat_HEAT2 = mat_HEAT
	mat_HEAT2[mat_HEAT2>1] <-1
	mat_HEAT2_ordered <- mat_HEAT[new_order,]
	heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>max(mat)*10],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20,(150*1/ncol(mat_HEAT2_ordered))),srtCol=45,cexRow=0.5,cexCol=0.5,main=i)

}

mat = t(table(droplevels(subset(ANVIO2, COG_FUNCTION  %in% c('Transposase', 'Integrase','transpos','integr'))[,c('combined','genome_name')])))
mat_HEAT = mat[,2:ncol(mat)]
	mat_HEAT2 = mat_HEAT
	mat_HEAT2[mat_HEAT2>1] <-1
	mat_HEAT2_ordered <- mat_HEAT[new_order,]
	heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>max(mat)*10],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20,(150*1/ncol(mat_HEAT2_ordered))),srtCol=45,cexRow=0.5,cexCol=0.5,main=i)


COG_list = c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W" )

mat = t(table(droplevels(subset(ANVIO2, COG_CATEGORY  == '' )[,c('protein_cluster_id','genome_name')])))
i="O"
	mat = t(table(droplevels(ANVIO2[grep(i, ANVIO2$COG_CATEGORY),c('combined','genome_name')])))
	mat_HEAT = mat[,2:ncol(mat)]
	mat_HEAT2 = mat_HEAT
	mat_HEAT2[mat_HEAT2>1] <-1
	mat_HEAT2_ordered <- mat_HEAT[new_order,]
	heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>10],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20,20),srtCol=45,cexRow=0.4,cexCol=0.4,main=i)


	heatmap.2(mat_HEAT2_ordered[,colSums(mat_HEAT2_ordered)>max(mat)*10],trace="none",col = inferno(75),Rowv= rep_tree_d,margins = c(20,(150*1/ncol(mat_HEAT2_ordered))),srtCol=45,cexRow=0.5,cexCol=0.5,main=i)


mat = t(table(droplevels(ANVIO2[grep('nitr', ANVIO2 $COG_FUNCTION),c('combined','genome_name')])))


selected = ANVIO2[grep('glutathione', ANVIO2$COG_FUNCTION),c('combined2','aa_sequence')]
D <- do.call(rbind, lapply(seq(nrow(selected)), function(i) t(selected[i, ])))


	write.table(D, file = '~/glutathione.fasta',row.names = FALSE, col.names = FALSE, quote = FALSE)

	write.table(D, file = paste('~/', phylogroup , '_CORE.fasta',sep=''),row.names = FALSE, col.names = FALSE, quote = FALSE)





JCORE = unique(COREmetadata[grep('J', COREmetadata$COG_CATEGORY),c("COG_FUNCTION")])
concatenated = c()
List_table = c()
for(i in JCORE){
	selected = ANVIO2[grep(i, 	ANVIO2$COG_FUNCTION),c('combined3','aa_sequence')]
	D <- do.call(rbind, lapply(seq(nrow(selected)), function(i) t(selected[i, ])))
	write.table(D, file = paste('~/', i , '_Q.fasta',sep=''),row.names = FALSE, col.names = FALSE, quote = FALSE)
	algn  = read.alignment(paste('~/', i , '_Q.fasta',sep=''), format = "fasta")
	if(length(unique(algn$nam))==60 && length(algn$nam) == 60){
		concatenated = cbind(concatenated,algn)
		algn2 = as.phyDat(algn,type="AA")
		dm  <- dist.ml(algn2)
		treeUPGMA  <- upgma(dm)
		treeNJ  <- NJ(dm)
		treeNJRooted =  midpoint.root(treeNJ)
		#increasing branch lengths
		treeNJRooted_sorted <- ladderize(treeNJRooted, right = FALSE)
 		plot(treeNJRooted_sorted, main=i)
		#plot(treeNJRooted_sorted, 'unrooted', main=i)
		write.table(D, file = paste('~/J/', i , '_Q.fasta',sep=''),row.names = FALSE, col.names = FALSE, quote = FALSE)
		List_table = rbind(List_table ,c(as.character(unique(ANVIO2[grep(i,ANVIO2$COG_FUNCTION),c('COG_FUNCTION_ACC')])[1]), i,length(attributes(algn2)$index),attributes(algn2)$nr))
	}
}

List_table = data.frame(List_table)
names(List_table) = c('COG','Function','aln_length','informativeSites')
write.csv(List_table, file = '~/Concatenated_table.csv')
List_table$aln_length = as.numeric(as.character( List_table$aln_length))
List_table$informativeSites = as.numeric(as.character( List_table$informativeSites))
plot(List_table$aln_length, List_table$informativeSites)


DCORE = unique(COREmetadata[grep('D', COREmetadata$COG_CATEGORY),c("COG_FUNCTION")])
concatenated = c()
List_table = c()
for(i in DCORE){
	selected = ANVIO2[grep(i, 	ANVIO2$COG_FUNCTION),c('combined3','aa_sequence')]
	D <- do.call(rbind, lapply(seq(nrow(selected)), function(i) t(selected[i, ])))
	write.table(D, file = paste('~/', i , '_Q.fasta',sep=''),row.names = FALSE, col.names = FALSE, quote = FALSE)
	algn  = read.alignment(paste('~/', i , '_Q.fasta',sep=''), format = "fasta")
	if(length(unique(algn$nam))==60 && length(algn$nam) == 60){
		concatenated = cbind(concatenated,algn)
		algn2 = as.phyDat(algn,type="AA")
		dm  <- dist.ml(algn2)
		treeUPGMA  <- upgma(dm)
		treeNJ  <- NJ(dm)
		treeNJRooted =  midpoint.root(treeNJ)
		#increasing branch lengths
		treeNJRooted_sorted <- ladderize(treeNJRooted, right = FALSE)
 		plot(treeNJRooted_sorted, main=i)
		#plot(treeNJRooted_sorted, 'unrooted', main=i)
		write.table(D, file = paste('~/D/', i , '_Q.fasta',sep=''),row.names = FALSE, col.names = FALSE, quote = FALSE)
		List_table = rbind(List_table ,c(as.character(unique(ANVIO2[grep(i,ANVIO2$COG_FUNCTION),c('COG_FUNCTION_ACC')])[1]), i,length(attributes(algn2)$index),attributes(algn2)$nr))
	}
}

List_table = data.frame(List_table)
names(List_table) = c('COG','Function','aln_length','informativeSites')
write.csv(List_table, file = '~/D_Concatenated_table.csv')
List_table$aln_length = as.numeric(as.character( List_table$aln_length))
List_table$informativeSites = as.numeric(as.character( List_table$informativeSites))
plot(List_table$aln_length, List_table$informativeSites)




Super_table = c()

COG_list = c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W" )
for(o in COG_list){
	DCORE = unique(COREmetadata[grep(o, COREmetadata$COG_CATEGORY),c("COG_FUNCTION")])
	concatenated = c()
	List_table = c()
	for(i in DCORE){
		selected = ANVIO2[grep(i, 	ANVIO2$COG_FUNCTION),c('combined3','aa_sequence')]
		D <- do.call(rbind, lapply(seq(nrow(selected)), function(i) t(selected[i, ])))
		write.table(D, file = paste('~/', i , '_Q.fasta',sep=''),row.names = FALSE, col.names = FALSE, quote = FALSE)
		algn  = read.alignment(paste('~/', i , '_Q.fasta',sep=''), format = "fasta")
		if(length(unique(algn$nam))==60 && length(algn$nam) == 60){
			concatenated = cbind(concatenated,algn)
			algn2 = as.phyDat(algn,type="AA")
			dm  <- dist.ml(algn2)
			treeUPGMA  <- upgma(dm)
			treeNJ  <- NJ(dm)
			treeNJRooted =  midpoint.root(treeNJ)
			#increasing branch lengths
			treeNJRooted_sorted <- ladderize(treeNJRooted, right = FALSE)
 			plot(treeNJRooted_sorted, main=i)
			#plot(treeNJRooted_sorted, 'unrooted', main=i)
			dir.create(o)
			write.table(D, file = paste('~/',o, '/', i , '_Q.fasta',sep=''),row.names = FALSE, col.names = FALSE, quote = FALSE)
			List_table = rbind(List_table ,c(as.character(unique(ANVIO2[grep(i,ANVIO2$COG_FUNCTION),c('COG_FUNCTION_ACC')])[1]), i,length(attributes(algn2)$index),attributes(algn2)$nr))}}
	List_table = data.frame(List_table)
	names(List_table) = c('COG','Function','aln_length','informativeSites')
	write.csv(List_table, file = paste('~/',o,'_Concatenated_table.csv',sep=''))
	List_table$aln_length = as.numeric(as.character( List_table$aln_length))
	List_table$informativeSites = as.numeric(as.character( List_table$informativeSites))
	plot(List_table$aln_length, List_table$informativeSites)
	Super_table =rbind(Super_table,List_table)
}
names(Super_table) = c('COG','Function','aln_length','informativeSites')
write.csv(Super_table, file = '~/Supertable.csv')
Super_table$aln_length = as.numeric(as.character( Super_table $aln_length))
Super_table$informativeSites = as.numeric(as.character( Super_table $informativeSites))
plot(Super_table$aln_length, Super_table $informativeSites)


install.packages('rlist')
library(rlist)


par(mfrow=c(2,5))
files = list.files(path = "~/COG_concatenation/")
listTrees = list()
plotlist = list()
for(file in files){
	algn  = read.alignment(paste('~/COG_concatenation/',file,sep=''), format = "fasta")
	concatenated = cbind(concatenated,algn)
	algn2 = as.phyDat(algn,type="AA")
	dm  <- dist.ml(algn2)
	treeUPGMA  <- upgma(dm)
	treeNJ  <- NJ(dm)
	treeNJRooted =  midpoint.root(treeNJ)
	#increasing branch lengths
	treeNJRooted_sorted <- ladderize(treeNJRooted, right = FALSE)
	plot(treeNJRooted_sorted, main=file)
	listTrees = list.append(listTrees, as.phylo(treeNJRooted_sorted))
	tr2 <- phylo4d(treeNJRooted_sorted, testGrouping,match.data = TRUE)
	t = ggtree(tr2) + geom_tippoint(aes( color= group), alpha=1)
	plotlist = list.append(plotlist,t)
	
}

consensusTree = ape::consensus(listTrees,p=0.5,check.labels=TRUE)
consensusTreeRooted_sorted <- ladderize(consensusTree, right = FALSE)

plot(consensusTreeRooted_sorted, main='Consensus')

plot(consensusTreeRooted_sorted, 'unrooted', main='Consensus')

algn  = read.alignment(paste('~/', i , '_Q.fasta',sep=''), format = "fasta")


testGrouping = ANVIO_cat[,2:3]
rownames(testGrouping) = c(as.character(ANVIO_cat[,1]))
require(phylobase)
tr2 <- phylo4d(consensusTreeRooted_sorted, testGrouping,match.data = TRUE)
ggtree(tr2) + geom_tiplab(aes(color=group)) +
    geom_tippoint(aes(shape= group, color= group), alpha=0.95)

consensusTree $edge.length[which(consensusTree $edge.length == 0)] <- 0.00001

tr2 <- phylo4d(consensusTree, testGrouping,match.data = TRUE)

ggtree(tr2,layout="fan") + geom_tippoint(aes( color= group), alpha=0.95)+theme(legend.position="right")
ggtree(tr2,layout="unrooted") + geom_tippoint(aes( color= group), alpha=0.95)+theme(legend.position="right")

par(mfrow=c(2,5))
t = ggtree(tr2) + geom_tippoint(aes( color= group), alpha=0.95)+theme(legend.position="right")

layout <-  matrix(c(1:25), nrow = 5, byrow = TRUE)
plotlist = list(t,t)
multiplot(plotlist = plotlist)


tr2 <- phylo4d(treeNJ, testGrouping,match.data = TRUE)

t3<-consensus.edges(listTrees,method="least.squares")


algn  = read.alignment('~/J/concatenation.fasta', format = "fasta")
algn2 = as.phyDat(algn)
dm  <- dist.ml(algn2)
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)
treeNJRooted =  midpoint.root(treeNJ)
#increasing branch lengths
treeNJRooted_sorted <- ladderize(treeNJRooted, right = FALSE)
plot(treeNJRooted_sorted, main=i)
plot(treeNJRooted_sorted,"unrooted", main="concatenated")




selected = ANVIO2[grep('DNA_gyrase-topoisomerase_IV-_subunit_B-DNA_gyrase-topoisomerase_IV-_subunit_B', ANVIO2$COG_FUNCTION),c('combined4','aa_sequence')]
D <- do.call(rbind, lapply(seq(nrow(selected)), function(i) t(selected[i, ])))
write.table(D, file = '~/GyrB.fasta',row.names = FALSE, col.names = FALSE, quote = FALSE)


algn  = read.alignment('~/Downloads/combined_aligned.fasta', format = "fasta")
algn = as.phyDat(algn)
dm  <- dist.ml(algn)
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)

 plot(treeNJ, "unrooted", main="NJ")
  plot(treeNJ, "unrooted", main="NJ")

  plot(treeUPGMA, "unrooted", main="NJ",cex=0.3)

 plot(treeNJ, main="NJ")
 plot(treeNJ, main="NJ")

treeNJRooted =  midpoint.root(treeNJ)
#increasing branch lengths
treeNJRooted_sorted <- ladderize(treeNJRooted, right = FALSE)
 plot(treeNJRooted_sorted, main="NJ")

 plot(treeNJRooted_sorted, 'unrooted', main="NJ")






ANVIO2$combined4 = paste(ANVIO2 $combined3, ANVIO2 $unique_id,sep="_ID_")

JCORE = unique(COREmetadata[grep('FoF1', COREmetadata$COG_FUNCTION),c("COG_FUNCTION")])
concatenated = c()
List_table = c()
for(i in JCORE){
	selected = ANVIO2[grep(i, 	ANVIO2$COG_FUNCTION),c('combined4','aa_sequence')]
	D <- do.call(rbind, lapply(seq(nrow(selected)), function(i) t(selected[i, ])))
	write.table(D, file = paste('~/', i , '_Q.fasta',sep=''),row.names = FALSE, col.names = FALSE, quote = FALSE)
	algn  = read.alignment(paste('~/', i , '_Q.fasta',sep=''), format = "fasta")
	concatenated = cbind(concatenated,algn)
	algn2 = as.phyDat(algn,type="AA")
	dm  <- dist.ml(algn2)
	treeUPGMA  <- upgma(dm)
	treeNJ  <- NJ(dm)
	treeNJRooted =  midpoint.root(treeNJ)
	#increasing branch lengths
	treeNJRooted_sorted <- ladderize(treeNJRooted, right = FALSE)
 	plot(treeNJRooted_sorted, main=i)
	#plot(treeNJRooted_sorted, 'unrooted', main=i)
	write.table(D, file = paste('~/FoF1/', i , '_Q.fasta',sep=''),row.names = FALSE, col.names = FALSE, quote = FALSE)
	List_table = rbind(List_table ,c(as.character(unique(ANVIO2[grep(i,ANVIO2$COG_FUNCTION),c('COG_FUNCTION_ACC')])[1]), i,length(attributes(algn2)$index),attributes(algn2)$nr))
}




selected = ANVIO2[grep('NAD-P-H-nitrite_reductase-_large_subunit-Rubredoxin', ANVIO2$COG_FUNCTION),c('combined4','aa_sequence')]
D <- do.call(rbind, lapply(seq(nrow(selected)), function(i) t(selected[i, ])))
write.table(D, file = '~/FoF1_gamma_subunit.fasta',row.names = FALSE, col.names = FALSE, quote = FALSE)


algn  = read.alignment('~/FoF1_gamma_subunit.fasta', format = "fasta")
algn = as.phyDat(algn)
dm  <- dist.ml(algn)
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)

plot(treeNJ, "unrooted", main="NJ")
plot(treeNJ, "unrooted", main="NJ")

plot(treeUPGMA, "unrooted", main="NJ",cex=0.3)

plot(treeNJ, main="NJ")
plot(treeNJ, main="NJ")

treeNJRooted =  midpoint.root(treeNJ)
#increasing branch lengths
treeNJRooted_sorted <- ladderize(treeNJRooted, right = FALSE)
plot(treeNJRooted_sorted, main="NJ")

plot(treeNJRooted_sorted, 'unrooted', main="NJ")



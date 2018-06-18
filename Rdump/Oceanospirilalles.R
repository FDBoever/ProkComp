##########################################################################
# Graphs and Statistics associated to CHECKM analysis
# FDB - March 2017 
# latest modifications - 29 April 2017 

######## PREPERATORY STEPS #########

# CheckM 3 should be run in nucleotide space (.fna) to accurately estimate the used parameters
# The -qa funtion of CheckM was converted in a matrix using Excell
# An additional column called "Group" (use environment or what ever type of grouping) was added, and assigned for each individual genome. 
# data matrix called qa_FNA.txt should be stored in your working directory
##########################################################################

######## LOAD AND INSTALL Libraries #########
library(ggplot2)
library(gridExtra)
library(rgl)
library(RColorBrewer)
library(Hmisc)
if(!require(dplyr)){install.packages("dplyr")}
if(!require(FSA)){install.packages("FSA")}
if(!require(DescTools)){install.packages("DescTools")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(multcompView)){install.packages("multcompView")}
library(igraph)
library(micropan)
library(reshape)
library(ade4)
library(dendextend)
library(igraph)
library(gplots)
library(RColorBrewer)
library(FactoMineR)
library(ggplot2)
library(Hmisc)


#########################
#
#
#########################

BactNet = readLines("~/DATA/MarinobacterGenomics/miscl/Bacterio_net_classifciation 2.txt")

cphylum = c()
cclass=c()
cOrder = c()
cFamily = c()
cgenus = c()

for(line in BactNet){
	list = strsplit(line, split=" ")
	if(list[[1]][1]=="Phylum"){
		phylum = list[[1]][2]
		}else if (list[[1]][1]=="Class"){
			class = list[[1]][2]
		}else if (list[[1]][1]=="Order"){
			order = list[[1]][2]
		}else if (list[[1]][1]=="Family"){
			family = list[[1]][2]
		}else{
			genus=list[[1]]
			genus = data.frame(genus)
			print(genus[1])
			cphylum = c(cphylum,rep(phylum,nrow(genus)))
			cclass = c(cclass,rep(class,nrow(genus)))
			cOrder = c(cOrder,rep(order,nrow(genus)))
			cFamily = c(cFamily,rep(family,nrow(genus)))
			cgenus = c(cgenus, as.character(genus$genus))
	}
}
TaxTable = cbind(Phylum = cphylum,Class= cclass,Order= cOrder,Family= cFamily,Genus_name= cgenus)
TaxTable=data.frame(TaxTable)

write.table(TaxTable,"TaxTable.txt",sep='\t')


########################
# 		Pfam 
########################


PfamData_out <- read.table('~/DATA/MarinobacterGenomics/miscl/PfamScan_output_oceano.csv',header=TRUE)
PfamData = table(PfamData_out$genome,PfamData_out$domain)

PfamData = PfamData[-grep("_HI0", rownames(PfamData)),]
PfamData = PfamData[-grep("Marinicella_sp", rownames(PfamData)),]

#count number of genomes in the dataset
number_of_genomes = length(colnames(PfamData))

# number of genomes per PC and Pfam
sumsPfam = data.frame(colSums(PfamData))
colnames(sumsPfam) = c('value')

#Combine in one dataframe

#Individual plots

ggplot(data= sumsPfam, aes(sumsPfam $value)) + 
  geom_histogram(bins = 60) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())


TotalDomains = rowSums(PfamData)
Core_Pfam_count = rowSums(PfamData[,colSums(PfamData)==60 ])
Pfam_singletons = PfamData[,colSums(PfamData)==1 ]
SummaryPfamPanGenome = data.frame(rowSums(Pfam_singletons))
colnames(SummaryPfamPanGenome) = c('Unique')
SummaryPfamPanGenome$name = rownames(SummaryPfamPanGenome)
SummaryPfamPanGenome$Core = Core_Pfam_count 
SummaryPfamPanGenome$Accessory = rowSums(PfamData) - SummaryPfamPanGenome$Core - SummaryPfamPanGenome$Unique

ggplot(data= SummaryPfamPanGenome, aes(y= SummaryPfamPanGenome$Unique,x= SummaryPfamPanGenome$name)) + geom_col() + xlab("Genome") + ylab("number of unique Pfam domains") + coord_flip() + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())


DF2 <- melt(SummaryPfamPanGenome, id.var="name")

ggplot(DF2, aes(x = name, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + xlab("Genome") + ylab("Pfam domains") + coord_flip() + scale_y_continuous(expand = c(0, 0)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank())



#----------------------------------------------
#	PAN GENOME CURVE
#----------------------------------------------
#---- species accumulation curve --------------

sp <- specaccum(PfamData, 'random', permutations=100)
Sp_Marinobacter <- specaccum( PfamData[grep("Marinobacter_",rownames(PfamData)),], 'random', permutations=100)
Sp_Alcanivorax <- specaccum( PfamData[grep("Alcanivorax",rownames(PfamData)),], 'random', permutations=100)
Sp_Marinobacterium <- specaccum( PfamData[grep("Marinobacterium",rownames(PfamData)),], 'random', permutations=100)
Sp_Halomonas <- specaccum( PfamData[grep("Halomonas",rownames(PfamData)),], 'random', permutations=100)
Sp_Marinomonas <- specaccum( PfamData[grep("Marinomonas",rownames(PfamData)),], 'random', permutations=100)



summary(sp)
plot(sp, ci.type='poly', col='darkred', lwd=2, ci.lty=0, ci.col='darkred', xlab='Genomes', ylab='Proteins', main='Gene accumulation plot')
boxplot(sp, col='white', add=TRUE,cex=0.5,pch = 21) 

plot(Sp_Halomonas, ci.type='poly', col='white', lwd=2, ci.lty=0, ci.col='white', xlab='Genomes', ylab='Domains', main='Gene accumulation plot',ylim=c(1500,4000))
#plot(sp, ci.type='poly', col='darkred', lwd=2, ci.lty=0, ci.col='darkred', xlab='Genomes', ylab='Proteins', main='Gene accumulation plot')

boxplot(Sp_Halomonas,  col = "white", border = "darkblue", add=TRUE,cex=0.5,pch = 21) 
boxplot(Sp_Marinobacter, col = "white", border = "orange", lty=1, cex=0.3, add= TRUE)
boxplot(Sp_Alcanivorax, col = "white", border = "black", lty=1, cex=0.3, add= TRUE)
boxplot(Sp_Marinobacterium, col = "white", border = "lightblue", lty=1, cex=0.3, add= TRUE)
boxplot(Sp_Marinomonas, col = "white", border = "red", lty=1, cex=0.3, add= TRUE)

plot(Sp_Halomonas, ci.type='poly', col='white', lwd=2, ci.lty=0, ci.col='white', xlab='Genomes', ylab='Domains', main='Gene accumulation plot',ylim=c(1500,4000))
#plot(sp, ci.type='poly', col='darkred', lwd=2, ci.lty=0, ci.col='darkred', xlab='Genomes', ylab='Proteins', main='Gene accumulation plot')

plot(Sp_Halomonas, ci.type='poly', col = "darkblue", add=TRUE,cex=0.3)
plot(Sp_Marinobacter,ci.type='poly',col = "orange", lty=2, cex=0.3, add= TRUE)
plot(Sp_Alcanivorax, ci.type='poly',col = "black", lty=2, cex=0.3, add= TRUE)
plot(Sp_Marinomonas, ci.type='poly',col = "red", lty=2, cex=0.3, add= TRUE)
plot(Sp_Marinobacterium, ci.type='poly',col = "lightblue", lty=2, cex=0.3, add= TRUE)


plot(Sp_Halomonas, ci.type='poly', col='black', lwd=2, ci.lty=0, ci.col='darkblue', xlab='Genomes', ylab='Proteins', main='accumulation plot',ylim=c(1750,4000))
plot(Sp_Marinobacter, ci.type='poly', col='black', lwd=2, ci.lty=0, ci.col='orange',add=TRUE)
plot(Sp_Alcanivorax, ci.type='poly', col='grey', lwd=2, ci.lty=0, ci.col='black',add=TRUE)
plot(Sp_Marinomonas, ci.type='poly', col='black', lwd=2, ci.lty=0, ci.col='red',add=TRUE)
plot(Sp_Marinobacterium, ci.type='poly', col='black', lwd=2, ci.lty=0, ci.col='lightblue',add=TRUE)



mods <- fitspecaccum(sp, "arrh")
modsHal <- fitspecaccum(Sp_Halomonas, "arrh")
modsMarb <- fitspecaccum(Sp_Marinobacter, "arrh")
modsAlc <- fitspecaccum(Sp_Alcanivorax, "arrh")
modsMarium <- fitspecaccum(Sp_Marinobacterium, "arrh")
modsMarinomon <- fitspecaccum(Sp_Marinomonas, "arrh")

plot(mods, col="grey",,xlab='Genomes', ylab='Protein Clusters', main='Gene accumulation plot')
plot(modsHal, col = "darkblue", add=TRUE,cex=0.3)
plot(modsMarb, col = "orange", add=TRUE,cex=0.3)
plot(modsAlc, col = "black", add=TRUE,cex=0.3)
plot(modsMarium, col = "lightblue", add=TRUE,cex=0.3)
plot(modsMarinomon, col = "red", add=TRUE,cex=0.3)



boxplot(sp, col = "white", border = "black", lty=1, cex=0.3, add= TRUE)
## Use nls() methods to the list of models
sapply(mods$models, AIC)

plot(sp, ci.type='poly', col='black', lwd=2, ci.lty=0, ci.col='lightskyblue3', xlab='Genomes', ylab='Proteins', main='accumulation plot')
boxplot(sp, col = "white", border = "darkblue", lty=1, cex=0.3, add= TRUE)

#---- sopen or closed panGENOME --------------
# Estimating if the pan-genome is open or closed based on a Heaps law model
#This function is based on a Heaps law approach suggested by Tettelin et al (2008). The Heaps law model is fitted to the number of new gene clusters observed when genomes are ordered in a random way. The model has two parameters, an intercept and a decay parameter called alpha. If alpha>1.0 the pan-genome is closed, if alpha<1.0< span=""> it is open.
#----------------------------------------------

heaps(PfamData,n.perm=100)
heaps(PfamData[grep("Marinobacter_",rownames(PfamData)),],n.perm=100)
heaps(PfamData[grep("Alcanivorax",rownames(PfamData)),],n.perm=100)
heaps(PfamData[grep("Marinobacterium",rownames(PfamData)),],n.perm=100)
heaps(PfamData[grep("Halomonas",rownames(PfamData)),],n.perm=100)
heaps(PfamData[grep("Marinomonas",rownames(PfamData)),],n.perm=100)

Sp_Marinobacter <- specaccum( PfamData[grep("Marinobacter_",rownames(PfamData)),], 'random', permutations=100)
Sp_Alcanivorax <- specaccum( PfamData[grep("Alcanivorax",rownames(PfamData)),], 'random', permutations=100)
Sp_Marinobacterium <- specaccum( PfamData[grep("Marinobacter_",rownames(PfamData)),], 'random', permutations=100)
Sp_Halomonas <- specaccum( PfamData[grep("Halomonas",rownames(PfamData)),], 'random', permutations=100)
Sp_Marinomonas <- specaccum( PfamData[grep("Marinomonas",rownames(PfamData)),], 'random', permutations=100)



#---- Pan-genome size --------------
#Chao - computes the Choa lower bound estimated number of gene clusters in a pan-genome

chao.Pfampansize <- chao(PfamData)
chao(PfamData[grep("Alcanivorax",rownames(PfamData)),])
chao(PfamData[grep("Halomonas",rownames(PfamData)),])
chao(PfamData[grep("Marinobacter_",rownames(PfamData)),])
chao(PfamData[grep("Marinobacterium",rownames(PfamData)),])
chao(PfamData[grep("Marinomonas",rownames(PfamData)),])

#---- Binomial mixture model --------------

binomix2 <- binomixEstimate(PfamData, K.range=2:15)
binomix3 <- binomixEstimate(PfamData[grep("Alcanivorax",rownames(PfamData)),], K.range=2:15)
binomix4 <- binomixEstimate(PfamData[grep("Halomonas",rownames(PfamData)),], K.range=2:15)
binomix5 <- binomixEstimate(PfamData[grep("Marinobacter_",rownames(PfamData)),], K.range=2:15)
binomix6 <- binomixEstimate(PfamData[grep("Marinobacterium",rownames(PfamData)),], K.range=2:15)
binomix7 <- binomixEstimate(PfamData[grep("Marinomonas",rownames(PfamData)),], K.range=2:15)

# Displaying the BIC.table

binomix2$BIC.table
# Summary of model

summary(binomix2)
summary(binomix3)
summary(binomix4)
summary(binomix5)
summary(binomix6)
summary(binomix7)

# count how many genomes used to calculate the
nrow(PfamData)
nrow(PfamData[grep("Alcanivorax",rownames(PfamData)),])
nrow(PfamData[grep("Halomonas",rownames(PfamData)),])
nrow(PfamData[grep("Marinobacter_",rownames(PfamData)),])
nrow(PfamData[grep("Marinobacterium",rownames(PfamData)),])
nrow(PfamData[grep("Marinomonas",rownames(PfamData)),])




#---- DISTANCES AND WEIGHTED DISTANCES --------------
#Manhattan and/or Jaccard distances between pan-genome profiles. Jaccard distance is based on similarity of shared genes, while Manhattan distance also includes similarity of lacking a certain gene (which is often recommended).(Snipen and Ussery, 2010)
#----------------------------------------------

#Jaccard distances based on panmatrix
Jdist.blast <- distJaccard(PfamData)

#manhatten distamces 
Mdist.blast <- distManhattan(PfamData)

#fluidityy
fluid.blast <- fluidity(PfamData)

#----------------------------------------------
#This weighting function computes weights for gene cluster according to their distribution in a pan-genome. 
#When computing distances between genomes or a PCA, it is possible to give weights to the different gene clusters, emphasizing certain aspects. 
#As proposed by Snipen & Ussery (2010), we have implemented two types of weighting: The default ‘"shell"’ type means gene families occuring frequently in the genomes, denoted shell-genes, are given large weight (close to 1) while those occurring rarely are given small weight (close to 0). The opposite is the ‘"cloud"’ type of weighting. Genes observed in a minority of the genomes are referred to as cloud-genes. Presumeably, the ‘"shell"’ weighting will give distances/PCA reflecting a more long-term evolution, since emphasis is put on genes who have just barely diverged away from the core. The ‘"cloud"’ weighting emphasizes those gene clusters seen rarely. Genomes with similar patterns among these genes may have common recent history. A ‘"cloud"’ weighting typically gives a more erratic or ‘noisy’ picture than the ‘"shell"’ weighting. 
#----------------------------------------------

PAPfamData = PfamData
PAPfamData[PAPfamData>1] <- 1 

w <- geneWeights(PfamData,type="shell")
w2 <- geneWeights(PfamData,type="cloud")
v <- geneWeights(PAPfamData,type="shell")
v2 <- geneWeights(PAPfamData,type="cloud")


Mdist.blast <- distManhattan(PfamData)
Jdist.blast <- distJaccard(PfamData)
Mdist1.blast <- distManhattan(PfamData,weights=w)
Mdist2.blast <- distManhattan(PfamData,weights=w2)

PAMdist.blast <- distManhattan(PAPfamData)
PAJdist.blast <- distJaccard(PAPfamData)
PAMdist1.blast <- distManhattan(PAPfamData,weights=w)
PAMdist2.blast <- distManhattan(PAPfamData,weights=w2)



m2 <- melt(as.matrix(Mdist.blast))[melt(upper.tri(as.matrix(Mdist.blast)))$value,]
m3 <- melt(as.matrix(Jdist.blast))[melt(upper.tri(as.matrix(Jdist.blast)))$value,]
m4 <- melt(as.matrix(Mdist1.blast))[melt(upper.tri(as.matrix(Mdist1.blast)))$value,]
m5 <- melt(as.matrix(Mdist2.blast))[melt(upper.tri(as.matrix(Mdist2.blast)))$value,]

m6 <- melt(as.matrix(PAMdist.blast))[melt(upper.tri(as.matrix(PAMdist.blast)))$value,]
m7 <- melt(as.matrix(PAJdist.blast))[melt(upper.tri(as.matrix(PAJdist.blast)))$value,]
m8 <- melt(as.matrix(PAMdist1.blast))[melt(upper.tri(as.matrix(PAMdist1.blast)))$value,]
m9 <- melt(as.matrix(PAMdist2.blast))[melt(upper.tri(as.matrix(PAMdist2.blast)))$value,]


constructedMatrix = data.frame(cbind(m2$value,m3$value,m4$value,m5$value,m6$value,m7$value,m8$value,m9$value))
colnames(constructedMatrix) = c("Manthattan_UW","Jaccard_UW","Manthattan_SW","Manthattan_CW","PA_Manthattan_UW","PA_Jaccard_UW","PA_Manthattan_SW","PA_Manthattan_CW")

q1 = ggplot(constructedMatrix, aes(Manthattan_UW, Jaccard_UW)) + geom_point() + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) + theme_classic()
q2 = ggplot(constructedMatrix, aes(Manthattan_UW, Manthattan_SW)) + geom_point() + theme_classic()
 #+ geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) 
q3 = ggplot(constructedMatrix, aes(Manthattan_UW, Manthattan_CW)) + geom_point() + theme_classic()
q4 = ggplot(constructedMatrix, aes(PA_Manthattan_UW, PA_Jaccard_UW)) + geom_point() + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) + theme_classic()
q5 = ggplot(constructedMatrix, aes(PA_Manthattan_UW, PA_Manthattan_SW)) + geom_point() + theme_classic()
 #+ geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) 
q6 = ggplot(constructedMatrix, aes(PA_Manthattan_UW, PA_Manthattan_CW)) + geom_point() + theme_classic()

q7 = ggplot(constructedMatrix, aes(Manthattan_UW,PA_Manthattan_UW)) + geom_point() + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) + theme_classic()
q8 = ggplot(constructedMatrix, aes(Manthattan_SW, PA_Manthattan_SW)) + geom_point() + theme_classic()
 #+ geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) 
q9 = ggplot(constructedMatrix, aes(Manthattan_CW, PA_Manthattan_CW)) + geom_point() + theme_classic()


multiplot(q1,q4,q7,q2,q5,q8,q3,q6,q9,cols=3)

# # Based on these results, we can reject the null hypothesis that these two matrices, spatial distance and ozone distance, are unrelated with alpha = .05. The observed correlation, r = 0.1636308, suggests that the matrix entries are positively associated.  So smaller differences in ozone are generally seen among pairs of stations that are close to each other than far from each other. Note that since this test is based on random permutations, the same code will always arrive at the same observed correlation but rarely the same p-value.
mantel.rtest(Mdist.blast, Jdist.blast, nrepet = 999)

DUFData = PfamData[,grep("DUF", colnames(PfamData))]

Pfamtree <- panTree(PfamData,nboot=100)
PAPfamtree <- panTree(PAPfamData,nboot=100)
DUFtree <- panTree(DUFData,nboot=100)


nnet <- neighborNet(Mdist.blast)


plot(as.phylo(Pfamtree$Htree),type="unrooted",no.margin = TRUE,cex=0.3)

xn <- as.network(Pfamtree$Htree) 
plot(as.phylo(hc), type = "unrooted", cex = 0.6,
     no.margin = TRUE)

par(mfrow=c(1,3)) 

plot(Pfamtree,xlab="Pfamtree",cex=0.3)
plot(PAPfamtree,xlab="PAPfamtree",cex=0.3)
plot(DUFtree,xlab="DUFPfamtree",cex=0.3)


w <- geneWeights(PfamData,type="shell")
v <- geneWeights(PAPfamData,type="shell")

wPfamtree = panTree(Pfamtree, weights=w,nboot = 100)
wPAPfamData = panTree(PAPfamData, weights=v,nboot = 100)

par(mfrow=c(1,2)) 
plot(tree2)
plot(wPfamtree,xlab="weighted Pfam Mantattan distance")
plot(wPAPfamData,xlab="weighted PAPfam domains Mantattan distance")
library(dendextend)
tanglegram(as.dendrogram(PAPfamtree$Htree),as.dendrogram(wPAPfamData $Htree))
tanglegram(as.dendrogram(PAPfamtree$Htree),as.dendrogram(DUFtree $Htree))

library(ape)
library(phytools)
Phylophlan = read.tree("~/DATA/MarinobacterGenomics/miscl/Oceano_clean.tree.nwk")
PhylophlanRooted =  midpoint.root(Phylophlan)
tree2 <- ladderize(PhylophlanRooted, right = FALSE)
tr = plot(tree2)
is_tip <- tree2$edge[,2] <= length(tree2$tip.label)
ordered_tips <- tree2$edge[is_tip, 2]
tree2$tip.label[ordered_tips]

library(DECIPHER)


force.ultrametric<-function(tree,method=c("nnls","extend")){
    method<-method[1]
    if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
        rooted=TRUE,trace=0)
    else if(method=="extend"){
        h<-diag(vcv(tree))
        d<-max(h)-h
        ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
            y=tree$edge[,2])
        tree$edge.length[ii]<-tree$edge.length[ii]+d
    } else 
        cat("method not recognized: returning input tree\n\n")
    tree
}

ult.nnls<-force.ultrametric(tree2)
tanglegram(ladderize(ult.nnls, right = FALSE), ladder(as.dendrogram(PAPfamtree$Htree), decreasing = TRUE))

PhyloDistMatrix<-cophenetic(tree2)
PhyloDistMatrix = PhyloDistMatrix[-grep("Marinicella_sp", rownames(PhyloDistMatrix)),-grep("Marinicella_sp", colnames(PhyloDistMatrix))]
PfamDist=dist(PfamData)
PfamDist  = PfamDist[ order(row.names(PfamDist)), ]
PfamDist  = PfamDist[ , order(colnames(PfamDist))]
DUFdist = distManhattan(DUFData)
PhyloDistMatrix  = PhyloDistMatrix[ order(row.names(PhyloDistMatrix)), ]
PhyloDistMatrix  = PhyloDistMatrix[ , order(colnames(PhyloDistMatrix))]
Mdist.blast = Mdist.blast[ , order(colnames(Mdist.blast))]
Jdist.blast = Jdist.blast[ , order(colnames(Jdist.blast))]
AAI  = AAI[ , order(colnames(AAI))]
AAI  = AAI[ order(row.names(AAI)), ]
OF  = OF[ , order(colnames(OF))]
OF  = OF[ order(row.names(OF)), ]
Og  = Og[ , order(colnames(Og))]
Og  = Og[ order(row.names(Og)), ]
#distTips(tree)
Mdist.blast
m2 <- melt(as.matrix(PfamDist))[melt(upper.tri(as.matrix(PfamDist)))$value,]
m3 <- melt(as.matrix(PhyloDistMatrix))[melt(upper.tri(as.matrix(PhyloDistMatrix)))$value,]
m4 <- melt(as.matrix(Mdist.blast))[melt(upper.tri(as.matrix(Mdist.blast)))$value,]
m5 <- melt(as.matrix(Jdist.blast))[melt(upper.tri(as.matrix(Jdist.blast)))$value,]
#m6 <- melt(as.matrix(PAJdist.blast))[melt(upper.tri(as.matrix(PAJdist.blast)))$value,]
m7 <- melt(as.matrix(DUFdist))[melt(upper.tri(as.matrix(DUFdist)))$value,]
m8=melt(as.matrix(AAI))[melt(upper.tri(as.matrix(AAI)))$value,]
m10 = melt(as.matrix(OF))[melt(upper.tri(as.matrix(OF)))$value,]
m11 = melt(as.matrix(Og))[melt(upper.tri(as.matrix(Og)))$value,]


DUFdist

constructedMatrix = data.frame(cbind(m2$value,m3$value,m4$value,m5$value,m7$value,m8$value,m10$value,m11$value))
colnames(constructedMatrix) = c("Pfam","phylogeny","PfamManhat",'JaccardPfam','DUFdist','AAI','Orthologous_fraction','nr_of_orthologs')

g1 = ggplot(constructedMatrix, aes(phylogeny, Pfam)) + geom_point(size=0.3) + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T)  +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()
g2 = ggplot(constructedMatrix, aes(phylogeny, PfamManhat)) + geom_point(size=0.3) + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T)  +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()
g3 = ggplot(constructedMatrix, aes(phylogeny, JaccardPfam)) + geom_point(size=0.3) + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T)  +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()
g5 = ggplot(constructedMatrix, aes(phylogeny, DUFdist)) + geom_point(size=0.3) + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()
g6 = ggplot(constructedMatrix, aes(phylogeny, 100-AAI)) + geom_point(size=0.3) + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()
multiplot(g6,g1,g2,g3,g5,cols=2)
g8 = ggplot(constructedMatrix, aes(PfamManhat, DUFdist)) + geom_point(size=0.3) + geom_smooth(method="lm", color="red", se=F, fullrange=T) +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()
g7 = ggplot(constructedMatrix, aes(phylogeny, 100-Orthologous_fraction)) + geom_point(size=0.3) + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()
g82 = ggplot(constructedMatrix, aes(phylogeny, nr_of_orthologs)) + geom_point(size=0.3) + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()
g9 = ggplot(constructedMatrix, aes(100-AAI, nr_of_orthologs)) + geom_point(size=0.3) + geom_smooth(method="lm", color="red", se=F, fullrange=T) +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()
g10 = ggplot(constructedMatrix, aes(AAI, Orthologous_fraction)) + geom_point(size=0.3) + geom_smooth(method="lm", color="red", se=F, fullrange=T) +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()
g11 = ggplot(constructedMatrix, aes(AAI, nr_of_orthologs)) + geom_point(size=0.3) + geom_smooth(method="lm", color="red", se=F, fullrange=T) +  geom_smooth(color="lightgrey", se=F, fullrange=T) + theme_classic()

mantel.rtest(PhyloDistMatrix, PfamDist, nrepet = 999)

multiplot(g6,g1,g3,g2,g5,g7,g8,g9, g82,g10,g11,cols=3)




testMatrix = cbind(m8$X1,m8$X2,constructedMatrix)
t2 = subset(testMatrix,phylogeny>2)
subset(t2, JaccardPfam <0.2)

pca3 = panpca(PfamData, scale = 0)
pca4 = panpca(PAPfamData, scale = 0)
pca5 = panpca(DUFData, scale = 0)
pca6 = panpca(AAI, scale = 0)
pca7 = panpca(data.frame(OF[,1:ncol(OF)-1]), scale = 0)
pca8 = panpca(Og, scale = 0)

Species=c()
for (name in rownames(PfamData)){
  print(name)
  print(regexpr("_", name)[1])
  print(substr(name,1, regexpr("_", name)[1]-1))
  Species = c(Species,substr(name,1, regexpr("_", name)[1]-1))
}

PfamCat = data.frame(cbind("genome"=rownames(PfamData),"Genus_name"=Species))
PfamCat <- join(PfamCat, TaxTable, by = "Genus_name")


Fam_cat = read.table("~/Pfam_oceano_sample.txt",header=TRUE)
names(Fam_cat) = c("sample_file_name","genome","family")

Genera=c()
for (genome in Fam_cat$genome){
  print(genome)
  print(regexpr("_", genome)[1])
  print(substr(genome,1, regexpr("_", genome)[1]-1))
  Genera = c(Genera,substr(genome,1, regexpr("_", genome)[1]-1))
}
Fam_cat$genera = Genera


matched = data.frame(genome= PfamCat $genome, family= Fam_cat[match(PfamCat $group, Fam_cat$genera), 3])

cols = colorRampPalette(c("red", "blue", "yellow","green","orange"))(length(unique(PfamCat$group)))
PfamCat$color <- factor(PfamCat$group, labels = cols)

PfamCat$fam <- matched$family
cols = colorRampPalette(c("red", "blue", "yellow","green","orange"))(length(unique(PfamCat$Family)))
PfamCat$color_fam <- factor(PfamCat$Family, labels = cols)


pcacolors <- PfamCat$color
pcacolors <- PfamCat$color_fam

names(pcacolors) <- rownames(PfamData)


plotScores(pca3, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
plotScores(pca4, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
plotScores(pca4, x = 2, y = 3, show.labels = FALSE, col = pcacolors, pch = 16)
plotScores(pca4, x = 3, y = 4, show.labels = FALSE, col = pcacolors, pch = 16)
plotScores(pca4, x = 4, y = 5, show.labels = FALSE, col = pcacolors, pch = 16)

plotScores(pca5, x = 1, y = 2, show.labels = FALSE, col = pcacolors, pch = 16)
plotScores(pca5, x = 2, y = 3, show.labels = FALSE, col = pcacolors, pch = 16)
plotScores(pca5, x = 3, y = 4, show.labels = FALSE, col = pcacolors, pch = 16)
plotScores(pca5, x = 4, y = 5, show.labels = FALSE, col = pcacolors, pch = 16)



Og = data.frame(cbind(Og, Species))
OF = data.frame(cbind(OF, Species))

Species=c()
for (name in rownames(Og)){
  print(name)
  print(regexpr("_", name)[1])
  print(substr(name,1, regexpr("_", name)[1]-1))
  Species = c(Species,substr(name,1, regexpr("_", name)[1]-1))
}
Og = data.frame(cbind(Og, as.factor(Species)))

Og.pca <- PCA(Og, quali.sup= ncol(Og),graph = FALSE)
par(mfrow=c(2,2))
plot(Og.pca , habillage = ncol(Og), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC1-2 Subspace",cex=2,centre = NULL)
plot(Pfm.pca, habillage = ncol(Pfm), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC2-3 Subspace",axes=2:3,cex=2,centre = NULL)
plot(Pfm.pca, habillage = ncol(Pfm), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC3-4 Subspace",axes=3:4,cex=2,centre = NULL)
plot(Pfm.pca, habillage = ncol(Pfm), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC4-5 Subspace",axes=4:5,cex=2,centre = NULL)

Og.pca <- PCA(Og)

~/Pfam_oceano_sample.txt

layout3d(matrix(1:6, 2, 3,byrow = TRUE),sharedMouse = TRUE)
plot3d(pca4$Scores[,1:3], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()
plot3d(pca4$Scores[,2:4], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()
plot3d(pca4$Scores[,3:5], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
plot3d(pca5$Scores[,1:3], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()
plot3d(pca5$Scores[,2:4], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()
plot3d(pca5$Scores[,3:5], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))

rgl.postscript("PCA_AAI_OF.pdf","pdf") 


legend3d("topright", legend = unique(PfamCat $fam), pch = 16, col = unique(PfamCat $color_fam), cex=1, inset=c(0.02))



layout3d(matrix(1:4, 2, 2,byrow = TRUE),sharedMouse = TRUE)
Sample.scaled<- data.frame(apply(AAI,2,scale))
Sample.scaled.2 <- data.frame(t(na.omit(t(Sample.scaled))))
pca.Sample.2 <- princomp(Sample.scaled.2)
plot3d(pca.Sample.2$scores[,1:3], col= pcacolors, size=12, type='t')
grid3d(c("x", "y+", "z"))
text3d(pca.Sample.2$Scores[,1:3],texts=rownames(Sample.scaled.2))
text3d(pca.Sample.2 $Scores[,1], pca.Sample.2 $Scores[,2], pca.Sample.2 $Scores[,3],texts=c(rownames(pca.Sample.2 $Scores)), cex= 0.7, pos=3)
next3d()
plot3d(pca.Sample.2$scores[,4:6], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()
Sample.scaled<- data.frame(apply(OF[,1:ncol(OF)-1],2,scale))
Sample.scaled.2 <- data.frame(t(na.omit(t(Sample.scaled))))
pca.Sample.2 <- princomp(Sample.scaled.2)
plot3d(pca.Sample.2$scores[,1:3], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()
plot3d(pca.Sample.2$scores[,4:6], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()
Sample.scaled<- data.frame(apply(Og[,1:ncol(Og)-1],2,scale))
Sample.scaled.2 <- data.frame(t(na.omit(t(Sample.scaled))))
pca.Sample.2 <- princomp(Sample.scaled.2)
plot3d(pca.Sample.2$scores[,1:3], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()
plot3d(pca.Sample.2$scores[,4:6], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()




plot3d(pca4$Scores[,1:3], col= pcacolors, size=20, type='p')
grid3d(c("x", "y+", "z"))
text3d(pca4$Scores[,1:3],texts=rownames(PfamData))


plot3d(pca3$Scores[,1:3], col= pcacolors, size=20, type='p')
grid3d(c("x", "y+", "z"))
legend3d("topright", legend = unique(PfamCat $group), pch = 16, col = unique(PfamCat $color), cex=1, inset=c(0.02))

layout3d(matrix(1:3, 1, 3),sharedMouse = TRUE)
plot3d(pca3$Scores[,1:3], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()
plot3d(pca3$Scores[,2:4], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))
next3d()
plot3d(pca3$Scores[,3:5], col= pcacolors, size=12, type='p')
grid3d(c("x", "y+", "z"))


plot3d(pca2$Scores[,1:3], col= pcacolors, size=20, type='p',main='test')
grid3d(c("x", "y+", "z"))
#text3d(pca2$Scores[,1], pca2$Scores[,2], pca2$Scores[,3],texts=c(rownames(pca2$Scores)), cex= 0.7, pos=3)
next3d()



legend("topright", col=unique(PfamCat $color), legend = unique(PfamCat $Genus_name),
    pch = 20, bty='n', cex=.75)

#---- BINOMIAL MIXTURE MODEL ------------------------------


STAMP1 = read.table(file = '~/DATA/MarinobacterGenomics/miscl/PA_vs_all_bonferroni_one_sided_CI_0.999_p_0.01.tsv', sep = '\t', header = TRUE)
hist(STAMP1$p.values..corrected.)
STAMP2 = subset(STAMP1,p.values..corrected.<0.05)
hist(STAMP2$p.values..corrected.)

STAMP2 = subset(STAMP1,p.values..corrected.<0.00001)
STAMP2 = STAMP1[grep("DUF", STAMP1$.),]
STAMP2 = subset(STAMP2,p.values..corrected.<0.00001)

STAMP2 = STAMP2[order(STAMP2$.), ]

hist(STAMP2$p.values..corrected.)
hist(STAMP2$Difference.between.means)

p2 = ggplot(STAMP2, aes(x = p.values..corrected., y = Difference.between.means)) + 
  geom_bar(stat = "identity",position="dodge") + coord_flip()
 
selectedSTAMP = data.frame("Mb_mean" = STAMP2[,2], "other_mean"= STAMP2[,4])#,"diff"=STAMP2[,8]))
rownames(selectedSTAMP)= STAMP2$.
selectedSTAMP$Pfam = STAMP2$.

melt(selectedSTAMP,id= selectedSTAMP$pfam)

Stamp_molten = melt(data.frame(cbind(as.character(STAMP2[,1]),STAMP2[,2],STAMP2[,4])),id=c("X1"))

as.numeric(paste(value))


ggplot(STAMP2, aes(x = p.values..corrected., y = Difference.between.means)) + 
  geom_bar(stat = "identity",position="dodge") + coord_flip()

grep(Stamp_molten$X1, "DUF")
Stamp_molten = subset(Stamp_molten, X1 ="DUF")
ggplot(Stamp_molten, aes(x = reorder(X1,value), y = as.numeric(paste(value)), fill = variable)) + 
  geom_bar(stat = "identity",position="dodge") + coord_flip()

ggplot(Stamp_molten, aes(x = reorder(X1,as.numeric(paste(value))), y = as.numeric(paste(value)), fill = variable)) + 
  geom_bar(stat = "identity",position="dodge") + coord_flip() + facet_wrap(~ X1)



ggplot(STAMP2, aes(x = p.values..corrected., y = Difference.between.means)) + 
  geom_bar(stat = "identity",position="dodge") + coord_flip()



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




nnet <- neighborNet(Mdist.blast)
#par("mar" = rep(1, 4))
plot(nnet, "2D")
plot(nnet)

########################
# SEQUENCE SIMMILARIY APPROAHCES
########################

TaxTable$genusX1 = TaxTable$Genus_name
TaxTable$genusX2 = TaxTable$Genus_name

#Make a corrected tax thingy, as we suspect that Marinobacter as well as Marinobacterium are Oceanospirillales
TaxTableCorrected = TaxTable
levels(TaxTableCorrected $Family)=c(levels(TaxTableCorrected $Family),"Marinobacteraceae",'-')
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacter",c("Order")]= "Oceanospirillales"
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacter",c("Family")]= "Marinobacteraceae"
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacterium",c("Order")]= "Oceanospirillales"
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacterium",c("Family")]= "-"


meltTax = function(corMat, metaTax){
	moltencorMat = melt(corMat)
	moltencorMat$species1 = gsub("^([^_]*_[^_]*)_.*$", "\\1", moltencorMat$X1)
	moltencorMat$species2 = gsub("^([^_]*_[^_]*)_.*$", "\\1", moltencorMat$X2)
	moltencorMat$genusX1 = gsub("\\_.*","", moltencorMat$X1)
	moltencorMat$genusX2 = gsub("\\_.*","", moltencorMat$X2)
	moltencorMat$genus_of_interest = ifelse(paste(moltencorMat$genusX1,"_",sep="")==GenusOfInterest,"yes","no")

	moltencorMat = merge((moltencorMat), metaTax[,c("genusX1", "Family", "Order")], by = 'genusX1')
	moltencorMat = merge((moltencorMat), metaTax[,c("genusX2", "Family", "Order")], by = 'genusX2')
	moltencorMat$X2 <- factor(moltencorMat$X2, levels=levels(moltencorMat$X1))
	moltencorMat$taxorder = ifelse(moltencorMat$X1==moltencorMat$X2,"itself",ifelse(moltencorMat$species1==moltencorMat$species2 & !(paste(moltencorMat$genusX1 ,"sp",sep="_") == moltencorMat$species1|paste(moltencorMat$genusX1 ,"sp.",sep="_") == moltencorMat$species1),"intraSpecies",ifelse(moltencorMat$genusX1==moltencorMat$genusX2,"intragenus",ifelse(moltencorMat$Family.x ==moltencorMat$Family.y,"intraFamily",ifelse(moltencorMat$Order.x == moltencorMat$Order.y,"intraOrder","intrerOrder")))))

	return(moltencorMat)
}




########################
# COMPARE M
########################
#install if needed
#install.packages("igraph",dependencies=TRUE)
#install.packages(c("igraph","gplots","RColorBrewer","FactoMineR","ggplot2","Hmisc"),dependencies=TRUE)


#load the CompareM data (Get rid of the spaces in headers in text-editor!)
CompareM =read.table("~/DATA/MarinobacterGenomics/miscl/AAI.tsv",header=TRUE)
#CompareM =read.table("~/DATA/MarinobacterGenomics/miscl/CompareM/AAI_oceanospirilalles.tsv",header=TRUE)

#show graphically
ggplot(data=CompareM,aes(x= Mean_AAI,y=Orthologous_fraction_OF))+geom_point()+geom_smooth(col="darkgreen",method="lm",se = TRUE)+theme_classic()
ggplot(data=CompareM,aes(x= Mean_AAI,y= X_orthologous_genes))+geom_point()+geom_point()+geom_smooth(col="darkgreen",method="lm",se = TRUE)+theme_classic()

# correlations
cor(CompareM$Mean_AAI, CompareM$orthologous_genes)
cor(CompareM$Mean_AAI, CompareM$Orthologous_fraction_.OF.)

# statistics of linear regression
summary(lm(CompareM$X_orthologous_genes ~ CompareM$Mean_AAI))
summary(lm(CompareM$Orthologous_fraction_OF ~ CompareM$Mean_AAI))

#transform the data in a matrix format
dat = CompareM[,c("Genome_A","Genome_B","Mean_AAI")]
g <- graph.data.frame(dat, directed=FALSE)
AAI = get.adjacency(g, attr="Mean_AAI", sparse=FALSE)

#genomes compared to themselves have 100 as a value
AAI[AAI<1] <- 100
AAI= AAI[-grep("Marinicella_", rownames(AAI)),-grep("Marinicella_", colnames(AAI))]



GenusOfInterest = "Marinobacter_"

moltenAAI = meltTax(AAI,TaxTable)
#moltenAAI$taxorder = as.factor(moltenAAI$taxorder)
#levels(moltenAAI$taxorder) = c("intrerOrder","intraOrder","intraFamily","intragenus","intraspecies","itself") 
g1.1 = ggplot(data= moltenAAI,aes(x= reorder(taxorder,value),y= value,colour= genus_of_interest))+geom_jitter(size=0.2)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10),axis.text.y = element_text(size=10))+scale_colour_manual(values=c("black","orange"))+geom_hline(yintercept = 95,colour='gray')+geom_hline(yintercept = 96,colour='gray')

g1.11 = ggplot(moltenAAI[moltenAAI $taxorder!='itself',], aes(x = value, fill = taxorder)) + geom_density(alpha = .5)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10),axis.text.y = element_text(size=10))+geom_vline(xintercept = 95,colour='gray')+geom_vline(xintercept = 96,colour='gray')+scale_fill_manual(values=c("black","orange",'red','yellow','green','blue'))+scale_y_continuous(expand = c(0, 0))

moltenAAI2 = meltTax(AAI, TaxTableCorrected)
g1.2 = ggplot(data= moltenAAI2,aes(x= reorder(taxorder,value),y= value,colour= genus_of_interest))+geom_jitter(size=0.2)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10),axis.text.y = element_text(size=10))+scale_colour_manual(values=c("black","orange"))+geom_hline(yintercept = 95,colour='gray')+geom_hline(yintercept = 96,colour='gray')

g1.21 = ggplot(moltenAAI2[moltenAAI2$taxorder!='itself',], aes(x = value, fill = taxorder)) + geom_density(alpha = .5)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10),axis.text.y = element_text(size=10))+geom_vline(xintercept = 95,colour='gray')+geom_vline(xintercept = 96,colour='gray')+scale_fill_manual(values=c("black","orange",'red','yellow','green','blue'))+scale_y_continuous(expand = c(0, 0))


multiplot(g1.1,g1.2,g1.11,g1.21,cols=2)
 
#heatmap
heatmap.2(AAI,trace="none",col = colorRampPalette(c("white","black"))(40))

#PCA
CompM.pca <- PCA(AAI)
plot(CompM.pca,cex=0.5)

#same, but for the Orhologous fraction
dat = CompareM[,c("Genome_A","Genome_B","Orthologous_fraction_OF")]
g <- graph.data.frame(dat, directed=FALSE)
OF = get.adjacency(g, attr="Orthologous_fraction_OF", sparse=FALSE)
OF[OF <1] <- 100
OF = OF[-grep("Marinicella_", rownames(OF)),-grep("Marinicella_", colnames(OF))]

heatmap.2(OF,trace="none",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))

CompM.pca <- PCA(OF)
plot(CompM.pca,cex=0.5)

#same, but for the Orhologous fraction
dat = CompareM[,c("Genome_A","Genome_B","X_orthologous_genes")]
g <- graph.data.frame(dat, directed=FALSE)
Og = get.adjacency(g, attr="X_orthologous_genes", sparse=FALSE)
Og[Og<1] <- NA
Og = Og[-grep("Marinicella_", rownames(Og)),-grep("Marinicella_", colnames(Og))]

heatmap.2(Og,trace="none",na.color="white",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))

CompM.pca <- PCA(Og)
plot(CompM.pca,cex=0.5)

AAI  = AAI[ , order(colnames(AAI))]
AAI  = AAI[ order(row.names(AAI)), ]
OF  = OF[ , order(colnames(OF))]
OF  = OF[ order(row.names(OF)), ]
Og  = Og[ , order(colnames(Og))]
Og  = Og[ order(row.names(Og)), ]

Og = data.frame(Og)
OF = data.frame(OF)
OF$Species=Species
Og$Species = Species

codon_usage.pca <- PCA(OF, quali.sup= ncol(OF))
par(mfrow=c(2,2))
plot(codon_usage.pca , habillage = ncol(OF), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC1-2 Subspace",cex=2)
plot(codon_usage.pca, habillage = ncol(codon_usage), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC2-3 Subspace",axes=2:3,cex=1)
plot(codon_usage.pca, habillage = ncol(Og), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC3-4 Subspace",axes=3:4,cex=1)


########################
# ANIb - pyani
########################
ANIb =read.table("~/DATA/MarinobacterGenomics/miscl/Genomes/pyani-out/ANIb_percentage_identity.tab",header=TRUE)
colnames(ANIb)=rownames(ANIb)
ANIb=ANIb[rownames(AAI),rownames(AAI)]

heatmap.2(as.matrix(ANIb),trace="none",col = colorRampPalette(c("white","black"))(40))
#ANIb = as.matrix(ANIb)

moltenANI = meltTax(as.matrix(ANIb), TaxTableCorrected)
moltenANI$taxorder = as.factor(moltenAAI$taxorder)
levels(moltenANI$taxorder) = c("intrerOrder","intraOrder","intraFamily","intragenus","intraspecies","itself") 
g1.1 = ggplot(data= moltenANI,aes(x= reorder(taxorder,value),y= value*100,colour= genus_of_interest))+geom_jitter(size=0.2)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10),axis.text.y = element_text(size=10))+scale_colour_manual(values=c("black","orange"))+geom_hline(yintercept = 95,colour='gray')+geom_hline(yintercept = 96,colour='gray')


moltenANI$AAI = moltenAAI2$value
ggplot(data= moltenANI,aes(x= value*100,y= AAI))+ geom_point(size=0.2)+xlab('ANIb')+theme_classic()+theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))

ggplot(data= moltenANI,aes(x= value*100,y= AAI))+ geom_point(size=0.2,aes(color=taxorder))+xlab('ANIb')+theme_classic()+theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+scale_colour_manual(values=colorRampPalette(c('#376268','#96a066'))(6))


ggplot(data= moltenANI,aes(x= value*100,y= AAI))+geom_point(size=0.2)+xlab('ANIb')+geom_smooth(method='lm',color='red')+geom_smooth()+theme_classic()+theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+facet_wrap(~taxorder,scales='free')






######### - statistics - ##########

# Mantel test 
#H0 - matrices are different, Ha - matrices correlate
mantel.rtest(as.dist(AAI), as.dist(OF), nrepet = 999)

orthologous_genes

#same, but for the Orhologous fraction

AAI = AAI[-grep("Marinicella_", rownames(AAI)),-grep("Marinicella_", colnames(AAI))]

Marinicella_sp._F2_strain_F2

AAI = AAI[grep("Marinobacter_", rownames(AAI)),grep("Marinobacter_", colnames(AAI))]
heatmap.2(AAI,trace="none",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))

aa_usage =read.table("~/DATA/MarinobacterGenomics/miscl/CompareM/aa_usage.txt",header=TRUE)
rownames(aa_usage)=aa_usage$Genome_ID
aa_usage = aa_usage[,2:ncol(aa_usage)]
heatmap.2(as.matrix(aa_usage),trace="none",na.color="white",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))

heatmap.2(as.matrix(aa_usage[grep("Marinobacter_", rownames(aa_usage)),colSums(aa_usage)>0.01]),trace="none",na.color="white",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40),margins=c(15,15))

aa_usage.pca <- PCA(aa_usage)
plot(aa_usage.pca,cex=0.5)


#####################################################################################################################

#####################################################################################################################
codon_usage =read.table("~/DATA/MarinobacterGenomics/miscl/CompareM/Oceanospir_codon_usage.txt",header=TRUE)
aa_usage =read.table("~/Downloads/aa_usage_1.txt",header=TRUE)

rownames(codon_usage)= codon_usage $Genome_ID
codon_usage = codon_usage[,2:ncol(codon_usage)]
codon_usage = codon_usage[-grep("_HI0", rownames(codon_usage)),]
codon_usage = codon_usage[,colSums(codon_usage)>100]


rownames(aa_usage)= aa_usage $Genome_ID
aa_usage = aa_usage[,2:ncol(aa_usage)]
aa_usage = aa_usage[-grep("_HI0", rownames(aa_usage)),]
#aa_usage = codon_usage[,colSums(aa_usage)>100]

m2 <- melt(as.matrix(vegdist(codon_usage)))[melt(upper.tri(as.matrix(vegdist(codon_usage))))$value,]
m3 <- melt(as.matrix(vegdist(aa_usage)))[melt(upper.tri(as.matrix(vegdist(aa_usage))))$value,]

constructedMatrix = data.frame(cbind(m2$value,m3$value))
colnames(constructedMatrix) = c("codon_usage","aa_usage")

ggplot(constructedMatrix, aes(codon_usage, aa_usage)) + geom_point() + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), color="red", se=F, fullrange=T) 


Species=c()
for (name in rownames(codon_usage)){
  print(name)
  print(regexpr("_", name)[1])
  print(substr(name,1, regexpr("_", name)[1]-1))
  Species = c(Species,substr(name,1, regexpr("_", name)[1]-1))
}

codon_usage$Species = Species
aa_usage$Species = Species

codom_usage.tmp = codon_usage
codon_usage.tmp$Species = NULL
heatmap.2(as.matrix(codon_usage.tmp),trace="none",na.color="white",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))#,RowSideColors= as.character(as.numeric(codon_usage$Species)))
heatmap.2(1-as.matrix(vegdist(codon_usage[,1:ncol(codon_usage)-1])),trace="none",na.color="white",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))#,RowSideColors=as.character(as.numeric(codon_usage$Species)))

codon_usage.pca <- PCA(codon_usage, quali.sup= ncol(codon_usage))
par(mfrow=c(2,2))
plot(codon_usage.pca , habillage = ncol(codon_usage), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC1-2 Subspace",cex=2)
plot(codon_usage.pca, habillage = ncol(codon_usage), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC2-3 Subspace",axes=2:3,cex=1)
plot(codon_usage.pca, habillage = ncol(codon_usage), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC3-4 Subspace",axes=3:4,cex=1)
plot(codon_usage.pca, habillage = ncol(codon_usage), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC4-5 Subspace",axes=4:5,cex=1)



codon_usage.pca <- PCA(Og, quali.sup= ncol(Og))
par(mfrow=c(2,2))
plot(codon_usage.pca , habillage = ncol(Og), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC1-2 Subspace",cex=2)
plot(codon_usage.pca, habillage = ncol(codon_usage), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC2-3 Subspace",axes=2:3,cex=1)
plot(codon_usage.pca, habillage = ncol(Og), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC3-4 Subspace",axes=3:4,cex=1)
plot(codon_usage.pca, habillage = ncol(Og), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC4-5 Subspace",axes=4:5,cex=1)




#####################################################################################################################
#####################################################################################################################

plot(as.igraph(codon_usage.pca, directed = TRUE, direction = "climb"))




codon_usage =read.table("~/DATA/MarinobacterGenomics/miscl/CompareM/codon_usage.txt",header=TRUE)
rownames(codon_usage)= codon_usage $Genome_ID
codon_usage = codon_usage[,2:ncol(codon_usage)]

codon_usage = codon_usage[-grep("_HI0", rownames(codon_usage)),]
codon_usage = codon_usage[,colSums(codon_usage)>100]

heatmap.2(as.matrix(codon_usage),trace="none",na.color="white",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))
heatmap.2(1-as.matrix(vegdist(codon_usage)),trace="none",na.color="white",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40), RowSideColors= as.character(as.numeric(codon_usage$Species)))


heatmap.2(as.matrix(codon_usage[grep("Marinobacter_", rownames(codon_usage)),colSums(codon_usage)>1]),trace="none",na.color="white",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))


codon_usage.pca <- PCA(codon_usage)
plot(codon_usage.pca,cex=0.5)

stop_usage =read.table("~/DATA/MarinobacterGenomics/miscl/CompareM/stop_usage.txt",header=TRUE)
rownames(stop_usage)= stop_usage $Genome_ID
stop_usage = stop_usage[,2:ncol(stop_usage)]
heatmap.2(as.matrix(stop_usage),trace="none",na.color="white",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))

stop_usage.pca <- PCA(stop_usage)
plot(stop_usage.pca,cex=0.5)



#######


source("http://bioconductor.org/biocLite.R") 
biocLite("FactoMineR")
biocLite("KARL")

#load packages
library(gplots)
library(FactoMineR)
library(vegan)
#load in 16S table
S16=read.table("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/16SSim_untrimmed.txt",header=TRUE)


#plot some heatmaps
heatmap.2(as.matrix(S16),trace="none")
heatmap.2(as.matrix(S16),trace="none",col = colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(40))
heatmap.2(as.matrix(S16),trace="none",col = colorRampPalette(c("red", "black", "yellow","darkgreen"))(40))
heatmap.2(as.matrix(S16),trace="none",col = colorRampPalette(c("white","grey", "black","yellow","darkgreen"))(40))

heatmap.2(as.matrix(S16),trace="none",col = colorRampPalette(c("#2c7fb8","#7fcdbb", "black","#edf8b1"))(40))


#PCA analysis with vegan
princomp(S16)
biplot(princomp(S16),display = "sites", xlabs=rep("●", nrow(S16)),cex=1)
 
 
#PCA analysis with factoMineR
#extract genus names from the organisms to later group them)
Species=c()
for (name in rownames(S16)){
  #print(name)
  #print(regexpr("_", name)[1])
  Species = c(Species,print(substr(name, regexpr("_", name)[1]+1, regexpr("_", name)[1]+4)))
  }
  
  
#add the genus names to the 16S data
S17 = cbind(S16, Species)
S16.pca <- PCA(S17, quali.sup= ncol(S17))
par(mfrow=c(2,2))
plot(S16.pca , habillage = ncol(S17), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC1-2 Subspace",cex=2)
plot(S16.pca, habillage = ncol(S17), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC2-3 Subspace",axes=2:3,cex=1)
plot(S16.pca, habillage = ncol(S17), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC3-4 Subspace",axes=3:4,cex=1)
plot(S16.pca, habillage = ncol(S17), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC4-5 Subspace",axes=4:5,cex=1)
#######

MAPLE_all_gamma = read.table("~/DATA/MarinobacterGenomics/miscl/MAPLE_all_gamma.txt",header=TRUE)
MAPLE_all_gamma_cat = read.table("~/DATA/MarinobacterGenomics/miscl/MAPLE_all_gamma_cat.txt",header=TRUE)


MAPLE_all_gammaS =  data.frame(cbind(colnames(MAPLE_all_gamma[,4:ncol(MAPLE_all_gamma)]), t(MAPLE_all_gamma[,4:ncol(MAPLE_all_gamma)])))
colnames(MAPLE_all_gammaS)=c("ID", MAPLE_all_gamma$ID)

JOINED <- right_join(MAPLE_all_gammaS, MAPLE_all_gamma_cat, by=c("ID"))

MAPLE <- data.frame(cbind(JOINED[,c(2:(ncol(JOINED)-5))],Species= JOINED$Genus))
MAPLE =na.omit(MAPLE)
Species=MAPLE$Species
Names = MAPLE
MAPLE[] <- lapply(MAPLE, function(x) as.numeric(as.character(x)))
MAPLE$Species=Species

MAPLE = MAPLE[!rownames(MAPLE) %in% c("786","779","780","800"), ]

dummy_marb=c()
for (i in MAPLE$Species){
	if(i == "Marinobacter"){
		dummy_marb = c(dummy_marb,"Marinobacter")
	}
	else{
		dummy_marb = c(dummy_marb,"Other") 
	}
}

MAPLE$Species = dummy_marb


MAPLE.pca <- PCA(MAPLE, quali.sup= 764)
ind.coord = MAPLE.pca$ind$coord


plot(MAPLE.pca, habillage = 764, 
     col.hab = c("black"))

plot(MAPLE.pca, habillage = 764, 
     col.hab = c("black"),axes=2:3)

plot(MAPLE.pca, habillage = 764, 
     col.hab = c("black"),axes=3:4)
     
plot(MAPLE.pca, habillage = 764, 
     col.hab = c("black"),axes=4:5)

MAPLE = MAPLE[!(MAPLE$Species %in% c("Carsonella","Portiera","Moranella","Evansia")),]

MAPLE.pca <- PCA(MAPLE, quali.sup= 764)
plot(MAPLE.pca, habillage = 764, 
     col.hab = c("black"))

plot(MAPLE.pca, habillage = 764, 
     col.hab = c("black"),axes=2:3)

plot(MAPLE.pca, habillage = 764, 
     col.hab = c("black"),axes=3:4)
     
plot(MAPLE.pca, habillage = 764, 
     col.hab = c("black"),axes=4:5)








MAPLE.pca <- PCA(MAPLE[!(MAPLE$Species %in% c("Carsonella","Portiera","Moranella","Evansia")),], quali.sup= 764)
plot(MAPLE.pca, habillage = 764, 
     col.hab = c("green", "blue", "red",'orange','purple'), 
     title = "Dataset projected onto PC1-2 Subspace")

MAPLE.pca <- PCA(MAPLE[(MAPLE$Species %in% c("Marinobacter","Alteromonas","Halomonas","Colwellia")),], quali.sup= 764)
plot(MAPLE.pca, habillage = 764, 
     col.hab = c("green", "blue", "red",'orange','purple'), 
     title = "Dataset projected onto PC1-2 Subspace")

MAPLE.pca <- PCA(MAPLE[(MAPLE$Species %in% c("Marinobacter","Alteromonas","Halomonas","Colwellia","Alcanivorax","Hahella","Pseudoalteromonas","Marinomonas","oleiphilus","Oleispira","Thalassolituus","Marinobacterium","Kushneria","Cobetia","Oceanospirillum","Marinospirillum","Neptunimonas","Bermanella","Oceanobacter","Kangiella","Nitricola","Idiomarina","Shewanella","Pseudomonas","Glaciecola","Agaribacter")),], quali.sup= 764)
plot(MAPLE.pca, habillage = 764, 
     col.hab = c("green", "blue", "red",'orange','purple'), 
     title = "Dataset projected onto PC1-2 Subspace",cex=2)
     
     
MAPLE2 = MAPLE[(MAPLE$Species %in% c("Marinobacter","Alteromonas","Halomonas","Colwellia","Alcanivorax","Hahella","Pseudoalteromonas","Marinomonas","oleiphilus","Oleispira","Thalassolituus","Marinobacterium","Kushneria","Cobetia","Oceanospirillum","Marinospirillum","Neptunimonas","Bermanella","Oceanobacter","Kangiella","Nitricola","Idiomarina","Shewanella","Pseudomonas","Glaciecola","Agaribacter")),]

MAPLE2=MAPLE2[,1:ncol(MAPLE2)-1]
MAPLE2= MAPLE2[, colSums(MAPLE2 != 0) > 0]

heatmap.2(as.matrix(MAPLE2[,1:ncol(MAPLE2)-1]),trace="none",RowSideColors=MAPLE$Species)


MAPLE= MAPLE[, colSums(MAPLE != 0) > 0]
MAPLE = MAPLE[(MAPLE$Species %in% c("Marinobacter","Alteromonas","Halomonas","Colwellia","Alcanivorax","Hahella","Pseudoalteromonas","Marinomonas","Oleiphilus","Oleispira","Thalassolituus","Marinobacterium","Kushneria","Cobetia","Oceanospirillum","Marinospirillum","Neptunimonas","Bermanella","Oceanobacter","Kangiella","Nitricola","Idiomarina","Shewanella","Glaciecola","Agaribacter")),]



heatmap.2(as.matrix(MAPLE[,1:ncol(MAPLE)-1]),trace="none")
heatmap.2(as.matrix(MAPLE[,1:ncol(MAPLE)-1]),trace="none",RowSideColors=as.character(as.numeric(MAPLE$Species)),col=colorRampPalette(c("black", "yellow"))(10))
legend("topright",      
    legend = unique(MAPLE$Species),
    col = unique(as.numeric(MAPLE$Species)), 
    lty= 1,             
    lwd = 5,           
    cex=.7
    )
    




MAPLE_all_gammaS =  MAPLE_all_gamma[,4:ncol(MAPLE_all_gamma)]
rownames(MAPLE_all_gammaS) = MAPLE_all_gamma$ID

MAPLE.pca <- PCA(t(MAPLE_all_gammaS))
plot(MAPLE.pca)

##########################################################################
#						PLOT GC 
##########################################################################
install.packages('zoo')

install.packages('seqinr')
library('seqinr')
library(zoo)


lna = read.fasta(file="~/Genomics/ProkComp/QC_Assemblies/Marinobacter_adhaerens_HP15_ASM16629.fna",as.string=TRUE)


library(zoo) # rollapply


lna <- read.fasta(file = files[1], seqtype = c("DNA"))

cat("# How many sequences\n")
length(lna)

# Access the n-th element of the list
n <- 1
summary(lna[n])

# sliding window
pdf(file="analysis/Rplots_sliding_window.pdf")
par(mfcol=c(1,1), cex=1.5, mai = c(1.2, 1.2, 0.1, 0.1)) # c(bottom, left, top, right)

windowsize <- 10000
x <- seq(from = 1, to = length(lna[[n]])-windowsize, by = windowsize) / 10^6
xlab <- "Position (Mbp)"

y <- rollapply(data = lna[[n]], width = windowsize, by = windowsize, FUN = GC)
plot(x, y, type="l", xlab=xlab, ylab="GC content")

GC.skew <- function(x){ y <- table(x); (y["c"] - y["g"]) / (y["c"] + y["g"]) }
y <- rollapply(data = lna[[n]], width = windowsize, by = windowsize, FUN = GC.skew)
plot(x, y, type="l", xlab=xlab, ylab="GC skew"); abline(h = 0)
plot(x, cumsum(y), type="l", xlab=xlab, ylab="Cumulative GC skew")

dev.off()

# Print R version and packages
sessionInfo()
##########################################################################
#						Pfam testing
##########################################################################
Pfm = read.table("~/Pfam_bacterioplankton_R.txt",header=TRUE)
Pfm = data.frame(t(Pfm))

pfpca = panpca(Pfm, scale = 0)
plot3d(pfpca$Scores[,1:3], size=20, type='p')
grid3d(c("x", "y+", "z"))

Pfm = Pfm[,grep("DUF", colnames(Pfm))]


Species=c()
for (name in rownames(Pfm)){
  print(name)
  print(regexpr("_", name)[1])
  print(substr(name,1, regexpr("_", name)[1]-1))
  Species = c(Species,substr(name,1, regexpr("_", name)[1]-1))
  }


Pfm[] <- lapply(Pfm, function(x) as.numeric(as.character(x)))
Pfm[] <- lapply(Pfm, function(x) as.numeric(as.character(x)))

Pfm = data.frame(cbind(Pfm, Species))

Pfm.pca <- PCA(Pfm, quali.sup= ncol(Pfm),graph = FALSE)
par(mfrow=c(2,2))
plot(Pfm.pca , habillage = ncol(Pfm), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC1-2 Subspace",cex=2,centre = NULL)
plot(Pfm.pca, habillage = ncol(Pfm), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC2-3 Subspace",axes=2:3,cex=2,centre = NULL)
plot(Pfm.pca, habillage = ncol(Pfm), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC3-4 Subspace",axes=3:4,cex=2,centre = NULL)
plot(Pfm.pca, habillage = ncol(Pfm), col.hab = c("green", "blue", "red",'orange','purple'), title = "Dataset projected onto PC4-5 Subspace",axes=4:5,cex=2,centre = NULL)

par(mfrow=c(2,2))
plotellipses(Pfm.pca, cex=2,centre = NULL)
plotellipses(Pfm.pca,axes=2:3,cex=2,centre = NULL)
plotellipses(Pfm.pca,axes=3:4,cex=2,centre = NULL)
plotellipses(Pfm.pca,axes=4:5,cex=2,centre = NULL)

Pfm[Pfm>0] <-1

Pfm = Pfm[(Pfm $Species %in% c("Marinobacter","Alteromonas","Halomonas","Colwellia","Alcanivorax","Hahella","Pseudoalteromonas","Marinomonas","Oleiphilus","Oleispira","Thalassolituus","Marinobacterium","Kushneria","Cobetia","Oceanospirillum","Marinospirillum","Neptunimonas","Bermanella","Oceanobacter","Kangiella","Nitricola","Idiomarina","Shewanella","Glaciecola","Agaribacter")),]


##########################################################################
#						TESTING TRAITOR
##########################################################################

TRAITAR = read.table("~/Genomics/Traitar/output2/phenotype_prediction/predictions_single-votes_combined.txt",header=TRUE)

######## LOAD THE DATA #########
CheckM = read.table("~/DATA/MarinobacterGenomics/miscl/qa2-1.txt",header=TRUE)

#extract some informative colums
multiData = cbind(CheckM$Genome_size,CheckM$GC,CheckM$nr_predicted_genes)
cols <- colorRampPalette(c("black","blue", "yellow","red"))(40)[CheckM$Group]
rownames(multiData) = CheckM$Organism
colnames(multiData) = c("Genome_size","GC","Predicted_genes")

rownames(CheckM2) = CheckM2$Bin_Id
CheckM2[grep("HI0", rownames(CheckM2)), ]


CheckM2 = subset(CheckM, Completeness > 98)
CheckM2 = CheckM2[order(CheckM2$Bin_Id),]


p1 <- ggplot(CheckM2, aes(x=reorder(Group, Genome_size, FUN=median), y=Genome_size)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))
p2 <- ggplot(CheckM2, aes(x=reorder(Group, nr_predicted_genes, FUN=median), y=nr_predicted_genes)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))
#p3 <- ggplot(CheckM2, aes(x=reorder(Group, Coding_density, FUN=median), y=Coding_density)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))
p4 <- ggplot(CheckM2, aes(x=reorder(Group, GC, FUN=median), y=GC)) + xlab("Genera") + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))

multiplot(p1, p2, p4, cols=1)

plot(CheckM2$Genome_size,CheckM2$nr_predicted_genes)

ggplot(CheckM2, aes(x= Genome_size, y= nr_predicted_genes,col=Group))+geom_point(aes(size=GC))+geom_text(aes(label= Bin_Id),hjust=0, vjust=0)
ggplot(CheckM2, aes(x= Genome_size, y= nr_predicted_genes,col=Group))+geom_point(aes(size=Completeness))+geom_text(aes(label= Bin_Id),hjust=0, vjust=0)
ggplot(CheckM2, aes(x= Genome_size, y= nr_predicted_genes,col=Group))+geom_point(aes(size= Contamination))+geom_text(aes(label= Bin_Id),hjust=0, vjust=0)

ggplot(CheckM2, aes(x= Genome_size, y= nr_predicted_genes,col=Group))+geom_point(aes(size= Contamination))

ggplot(subset(CheckM2, Contamination<15), aes(x= Genome_size, y= nr_predicted_genes,col=Group))+geom_point(aes(size= Contamination))


CheckM2 = subset(CheckM2, Contamination<15)



p4 <- ggplot(CheckM2[grep("HI0", rownames(CheckM2)), ], aes(x=reorder(Group, GC, FUN=median), y=GC)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))


q = ggplot(CheckM2, aes(x=reorder(Group,GC), y=GC)) + 
    geom_boxplot() +      # Thinner lines
    xlab("Genome") +
    ylab("GC content") +
    theme_bw() + coord_flip()

q1 = ggplot(CheckM2, aes(x=Group, y=0.000001*Genome_size)) + 
    geom_boxplot() +      # Thinner lines
    xlab("Genome") +
    ylab("GC content") +
    theme_bw() + coord_flip()

q2 = ggplot(CheckM2, aes(x=Group, y= nr_predicted_genes)) + 
    geom_boxplot() +      # Thinner lines
    xlab("Genome") +
    ylab("GC content") +
    theme_bw() + coord_flip()




ggplotRegression <- function (fit) {

require(ggplot2)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "darkgreen") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))+theme_classic()
}

ggplotRegression(lm(nr_predicted_genes ~ Genome_size, data = CheckM2))




##########################################################################
#						COMPARE ANVIO AND CHECKM
##########################################################################
library(RColorBrewer)




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










##########################################################################
#						STATISTICAL ANALYSIS
##########################################################################


######## Pearson's Correlation Test #########

#Pearson correlation test
rcorr(multiData, type="pearson")


######## ANALYSIS OF VARIANCE #########

# ONE-WAY ANOVA

# testing assumptions
tapply(CheckM$GC, CheckM$Group, var)
#test for heteroscedasticity using the Fligner-Killeen test of homogeneity of variance
fligner.test(CheckM$GC~CheckM$Group)
#there was no evidence of any significant differences in variance accros samples, so it's legimitate to continue with one-way analysis of variance
#in case of GC (is significant, so look for alternative)
attach(CheckM)
summary(aov(GC~Group))
par(mfrow=c(2,2))
plot(aov(GC~Group))

#NON-parametric alternative (Kruskal-Wallis)
library(dunn.test)
kruskal.test(GC ~ Group, data = CheckM) 
PT = dunn.test(CheckM$GC, CheckM$Group,method="bh") #Benjamini-Hochberg correction
PT


##########################################################################
#						GRAPHS AND FIGURES
##########################################################################

######## VALUABLE FUNCTIONS #########
# multiplot() will allow you to visualize multiple plots in a grid arrangment

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

lm_eqn = function(m) {

  l <- list(a = format(coef(m)[1], digits = 2),
      b = format(abs(coef(m)[2]), digits = 2),
      r2 = format(summary(m)$r.squared, digits = 3));

  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }

  as.character(as.expression(eq));                 
}

######## boxplots with Group grouping #########





######## BAR PLOTS #########

q = ggplot(CheckM, aes(x=Organism, y=GC)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=.3) +      # Thinner lines
    geom_errorbar(aes(ymin= GC-GC_std, ymax= GC +GC_std),
                  size=.3,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("Genome") +
    ylab("GC content") +
    ggtitle("Summary statistics") +
    theme_bw()
q + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5)) + coord_cartesian(ylim = c(45, 70))



######## SCATTERPLOTS #########

# add's a linear model fit (geom_smooth) and confidence interval (se = TRUE) ()
p5 = ggplot(CheckM, aes(Genome_size,nr_predicted_genes,color = cols)) + geom_point() + geom_smooth(method="lm",color = "orange", size = 0.75) + geom_text(aes(x = 25, y = 300, label = lm_eqn(lm(nr_predicted_genes ~ Genome_size, CheckM))), parse = TRUE)

p6 = ggplot(CheckM, aes(Genome_size,nr_predicted_genes, col = Group)) + geom_point() + geom_smooth(method="lm", size = 0.75)
p7 = ggplot(CheckM, aes(Genome_size,nr_predicted_genes, col = Group)) + geom_point() + geom_smooth(method="lm", size = 0.75,se=FALSE,fullrange=TRUE)
p8 = ggplot(CheckM, aes(Genome_size,GC,color=Group)) + geom_point() + geom_smooth(method="lm",color = "orange", size = 0.75)
p9 = ggplot(CheckM, aes(Genome_size, GC, col = Group)) + geom_point() 
p10 = ggplot(CheckM, aes(Genome_size,Coding_density)) + geom_point() + geom_smooth(method="lm",color = "orange", size = 0.75)
p11 = ggplot(CheckM, aes(Genome_size, Coding_density, col = Group)) + geom_point() + geom_smooth(method="lm", size = 0.75,se=FALSE,fullrange=TRUE)
p12= ggplot(CheckM, aes(Genome_size, Contamination)) + geom_point() + geom_smooth(method="lm",color = "orange", size = 0.75)
p13= ggplot(CheckM, aes(Genome_size, Contamination, col = Group)) + geom_point() + geom_smooth(method="lm", size = 0.75,se=FALSE,fullrange=TRUE)

multiplot(p5, p6, p7, p8,p9,p10,p11,p12,p13, cols=3)
multiplot(p5, p8, cols=1)


######## ROTATABLE 3D PLOT #########

cols <- palette(brewer.pal(11,'Set3'))[CheckM$Group]

plot3d(multiData, col=cols,size=2,pch = 16)
rgl.postscript("3Dplot_MB_statistics.pdf","pdf")








######## IGRAPH #########
g1 <- graph( edges=c(1,2, 2,3, 3, 1), n=3, directed=F ) 
plot(g1) 
g2 <- graph( edges=c(1,2, 2,3, 3, 1), n=10 )

plot(g2)   
tr <- make_tree(40, children = 3, mode = "undirected")

plot(tr, vertex.size=10, vertex.label=NA) 

st <- make_star(40)

plot(st, vertex.size=10, vertex.label=NA) 

er <- sample_gnm(n=100, m=40) 

plot(er, vertex.size=6, vertex.label=NA)  


 ba <-  sample_pa(n=100, power=1, m=1,  directed=F)

 plot(ba, vertex.size=6, vertex.label=NA)

 zach <- graph("Zachary") # the Zachary carate club

 plot(zach, vertex.size=10, vertex.label=NA)

net <- graph_from_data_frame(d=codon_usage, vertices=nodes, directed=T) 


######## boxplots with Group grouping #########
links2 <- as.matrix(AAI)
net2 <- graph_from_incidence_matrix(links2)

# A built-in vertex attribute 'type' shows which mode vertices belong to.
table(V(net2)$type)

plot(net2,vertex.label=NA)

# To transform a one-mode network matrix into an igraph object,
# use graph_from_adjacency_matrix()

# We can also easily generate bipartite projections for the two-mode network:
# (co-memberships are easy to calculate by multiplying the network matrix by
# its transposed matrix, or using igraph's bipartite.projection function)

net2.bp <- bipartite.projection(net2)

# We can calculate the projections manually as well:
#   as_incidence_matrix(net2)  %*% t(as_incidence_matrix(net2))
# t(as_incidence_matrix(net2)) %*%   as_incidence_matrix(net2)

plot(net2.bp$proj1, vertex.label.color="black", vertex.label.dist=1,
     vertex.label=nodes2$media[!is.na(nodes2$media.type)])

plot(net2.bp$proj2, vertex.label.color="black", vertex.label.dist=1,
     vertex.label=nodes2$media[ is.na(nodes2$media.type)])
     
source("https://bioconductor.org/biocLite.R")
biocLite()
 
     
source("https://bioconductor.org/biocLite.R")
biocLite("FindMyFriends")

     
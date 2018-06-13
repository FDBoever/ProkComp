#MAPLE

MAPLE <- read.table('~/DATA/MarinobacterGenomics/2018_ProkComp/MAPLE_v1.txt',header=TRUE)

MAPLE2 = MAPLE[,c(6:length(colnames(MAPLE)))]
rownames(MAPLE2) = paste(MAPLE$Name,rownames(MAPLE),sep='_')


MAPLE3 = MAPLE2
MAPLE3$max = apply(MAPLE3, 1, max)
MAPLE3$min = apply(MAPLE3, 1, min)


MAPLE4 = MAPLE3[MAPLE3$max == 100,]
MAPLE4 = MAPLE4[MAPLE4$min < 100,c(1:(length(colnames(MAPLE4))-2))]


#Remove all modules which are absent in set of genomes

row_sub = apply(MAPLE2, 1, function(row) all(row !=0 ))
subset = MAPLE2[row_sub,]

heatmap.2(as.matrix(subset),trace='none',cexRow=0.3,cexCol=.5,col
=colorRampPalette(c('red', "yellow", "forestgreen"))(n = 299))

heatmap.2(as.matrix(subset),trace='none',cexRow=0.3,cexCol=.5,col
=colorRampPalette(c('white', "yellow", "forestgreen"))(n = 100))


MAPLEselect = MAPLE[MAPLE$Category =='Two-component_regulatory_system',]


MAPLE2 = MAPLEselect[,c(6:length(colnames(MAPLEselect)))]
rownames(MAPLE2) = paste( MAPLEselect $Name,rownames(MAPLEselect),sep='_')

row_sub = apply(MAPLE2, 1, function(row) all(row !=0 ))
subset = MAPLE2[row_sub,]

heatmap.2(as.matrix(MAPLE2),trace='none',cexRow=0.3,cexCol=.5,col
=colorRampPalette(c('red', "yellow", "forestgreen"))(n = 299))

heatmap.2(as.matrix(MAPLE2),trace='none',cexRow=0.5,cexCol=.5,col
=colorRampPalette(c('white', "yellow", "forestgreen"))(n = 100),margins=c(10,10))


#############################################################################

unique(MAPLE$Classification_A)
MAPLEselect = MAPLE[MAPLE$Category =='Nitrogen_metabolism',]
MAPLEselect = MAPLE[MAPLE$Category =='Nucleotide_and_amino_acid_metabolism',]
Phosphotransferase_system_(PTS) 

MAPLE2 = MAPLEselect[,c(6:length(colnames(MAPLEselect)))]
rownames(MAPLE2) = paste(MAPLEselect $Category, MAPLEselect $Name,rownames(MAPLEselect),sep='_')

row_sub = apply(MAPLE2, 1, function(row) all(row !=0 ))
subset = MAPLE2[row_sub,]

heatmap.2(as.matrix(subset),trace='none',cexRow=0.3,cexCol=.5,col
=colorRampPalette(c('red', "yellow", "forestgreen"))(n = 299))

heatmap.2(as.matrix(subset),trace='none',cexRow=0.5,cexCol=.5,col
=colorRampPalette(c('white', "yellow", "forestgreen"))(n = 100),margins=c(12,20))



heatmap.2(as.matrix(MAPLE2),trace='none',cexRow=0.5,cexCol=.5,col
=colorRampPalette(c('white', "yellow", "forestgreen"))(n = 100),margins=c(12,20))


for(category in unique(MAPLE$Category)){
	tryCatch({
	MAPLEselect = MAPLE[MAPLE$Category == category,]
	MAPLE2 = MAPLEselect[,c(6:length(colnames(MAPLEselect)))]
rownames(MAPLE2) = paste(MAPLEselect $Category, MAPLEselect $Name,rownames(MAPLEselect),sep='_')
	MAPLE2 = MAPLE2[,colnames(MAPLE2) %in%CheckMANVIO2$Bin_Id]
	MAPLE2 = MAPLE2[apply(MAPLE2, 1, max) ==100,]
	pdf(paste("maple_", category, ".pdf", sep = ""))
	heatmap.2(as.matrix(MAPLE2),trace='none',cexRow=0.5,cexCol=.5,col
=colorRampPalette(c('white', "yellow", "forestgreen"))(n = 100),margins=c(12,12))
	dev.off()
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


	MAPLE2 = MAPLE[,c(6:length(colnames(MAPLE)))]
	rownames(MAPLE2) = paste(MAPLEselect $Name,rownames(MAPLE),sep='_')

	MAPLE2 = MAPLE2[,colnames(MAPLE2) %in%CheckMANVIO2$Bin_Id]
	MAPLE2 = MAPLE2[apply(MAPLE2, 1, max) ==100,]
	heatmap.2(as.matrix(MAPLE2),trace='none',cexRow=0.5,cexCol=.5,col
=colorRampPalette(c('white', "yellow", "forestgreen"))(n = 100),margins=c(12,12))


Maple_tree = drop.tip(RAxMLCampbell, RAxMLCampbell $tip.label[-match(colnames(MAPLE2), RAxMLCampbell $tip.label)])

Maple_tree =  midpoint.root(Maple_tree)

#increasing branch lengths
Maple_tree <- ladderize(Maple_tree, right = FALSE)

#can’t have any branch lengths of zero or downstream commands will collapse those nodes…
Maple_tree $edge.length[which(Maple_tree $edge.length == 0)] <- 0.00001
Maple_tree_um <- chronopl(Maple_tree,lambda = 0.1,tol = 0)
Maple_tree_d <- as.dendrogram(as.hclust.phylo(Maple_tree_um))


#force row order so that it matches the order of leafs in rep_tree_d
clade_order <- order.dendrogram(Maple_tree_d)
clade_name <- labels(Maple_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, colnames(MAPLE2))
MAPLE2 <- MAPLE2[,new_order]






MAPLEselect = MAPLE[]
MAPLE2 = MAPLEselect[,c(6:length(colnames(MAPLEselect)))]
rownames(MAPLE2) = paste(MAPLEselect $Category, MAPLEselect $Name,rownames(MAPLEselect),sep='_')
MAPLE2 = MAPLE2[,colnames(MAPLE2) %in%CheckMANVIO2$Bin_Id]
MAPLE2 = MAPLE2[apply(MAPLE2, 1, max) ==100,]
MAPLE2 = MAPLE2[!rowSums(MAPLE2)/dim(MAPLE2)[2] ==100,]

heatmap.2(as.matrix(MAPLE2),trace='none',cexRow=0.5,cexCol=.5,col
=colorRampPalette(c('white', "yellow", "forestgreen"))(n = 100),margins=c(12,12))

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

heatmap.2(AAI,trace="none",scale="none",col=colorRampPalette(c("blue", 'black', "yellow", "red"))(n = 15))
heatmap.2(AAI,trace="none",scale="none",col=colorRampPalette(c("darkblue", 'teal', "yellow"))(n = 15))


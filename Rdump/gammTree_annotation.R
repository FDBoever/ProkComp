#
#	gammaTree_annotation.R
#
# 	Script designed to flexibly visualise large tree, example a Gammaproteobacterial tree with +1000 organisms
#
####################

library(ggtree)
library(dplyr)
library(phytools)
library(colorspace)
#####################


bigTree = read.tree('~/DATA/MarinobacterGenomics/miscl/FAA_Gamma_MB2.tree.nwk')
TaxTable = read.table('~/DATA/MarinobacterGenomics/miscl/TaxTable.txt')


#bigTree = read.tree('~/DATA/MarinobacterGenomics/miscl/rerooted_tree')

#Define outgroup
outgroup <- c("Phaeobacter_inhibens_2.10", "Phaeobacter_inhibens_DSM_17395","Phaeobacter_gallaeciensis_DSM_26640",
"Phaeobacter_sp._P97","Ruegeria_sp._TM1040","Ruegeria_mobilis","Ruegeria_pomeroyi","Roseobacter_litoralis","Roseobacter_denitrificans","Dinoroseobacter_shibae","Rhodobacter_sphaeroides_2.4.1","Rhodobacter_sphaeroides_KD131","Rhodobacter_sphaeroides_ATCC_17029","Rhodobacter_sphaeroides_ATCC_17025","Rhodobacter_sp._LPB0142","Rhodobacter_capsulatus","Pseudomonas_geniculata","Alpha_proteobacterium_HIMB5","Alpha_proteobacterium_HIMB59") # I choose this clade !

rooted_tree <- root(bigTree, outgroup, resolve.root=TRUE)

midRoot_bigTree = ladderize(rooted_tree)

genera_tip_labels = gsub("\\_.*","", midRoot_bigTree$tip.label)
genera_tip_labels = data.frame(genera_tip_labels)
colnames(genera_tip_labels)=c('Genus_name')



TaxTable$genusX1 = TaxTable$Genus_name
TaxTable$genusX2 = TaxTable$Genus_name

#Make a corrected tax thingy, as we suspect that Marinobacter as well as Marinobacterium are Oceanospirillales
TaxTableCorrected = TaxTable
levels(TaxTableCorrected $Family)=c(levels(TaxTableCorrected $Family),"Marinobacteraceae",'-')
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacter",c("Order")]= "Oceanospirillales"
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacter",c("Family")]= "Marinobacteraceae"
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacterium",c("Order")]= "Oceanospirillales"
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacterium",c("Family")]= "-"


taxa.ordered = right_join(x = TaxTableCorrected, y = genera_tip_labels, by = "Genus_name")

levels(taxa.ordered$Phylum)=c(levels(taxa.ordered$Phylum),'unknown')
levels(taxa.ordered$Class)=c(levels(taxa.ordered$Class),'unknown')
levels(taxa.ordered$Order)=c(levels(taxa.ordered$Order),'unknown')
levels(taxa.ordered$Family)=c(levels(taxa.ordered$Family),'unknown')
levels(taxa.ordered$Phylum)=c(levels(taxa.ordered$Phylum),'unknown')
levels(taxa.ordered$genusX1)=c(levels(taxa.ordered$genusX1),'unknown')
levels(taxa.ordered$genusX2)=c(levels(taxa.ordered$genusX2),'unknown')

taxa.ordered[is.na(taxa.ordered)] <- "unknown" 

rownames(taxa.ordered)=midRoot_bigTree$tip.label
taxa.ordered$strain = midRoot_bigTree$tip.label

l <- list(); 
for(Order in unique(taxa.ordered$Order)){
l[[Order]] = as.character(taxa.ordered[taxa.ordered $Order==Order,'strain'])
}

annotatedTree <- groupOTU(midRoot_bigTree, l)


ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,) +geom_treescale()+theme_tree() +  scale_color_manual(values=c("black", rainbow_hcl(17))) + theme(legend.position="right")


p6= ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,) +geom_treescale()+theme_tree() +  scale_color_manual(values=rev(colorRampPalette(brewer.pal(8, "Set1"))(18))) + theme(legend.position="right")

colint = 200
	for(genus in  c('Marinobacter','Alteromonas','Escherichia', 'Vibrio','Salmonella', 'Shewanella','Acinetobacter','Yersinia','Francisella','Buchnera','Alcanivorax','Legionella','Xanthomonas','Haemophilus')){
	genusTip = annotatedTree $tip.label[grepl(paste(genus,'_',sep=''), annotatedTree $tip.label)]
	rca = getMRCA(annotatedTree,tip= genusTip)
	print(rca)
	p6 = p6 + geom_hilight(node= rca, fill=colors()[colint])+geom_cladelabel(node= rca, label= genus,offset.text=0.3)
	colint = colint + 10 
}

p6+xlim(0,9)






############################################################

CAZYf = read.csv('~/DATA/MarinobacterGenomics/miscl/Genomes/dbCAN/dbCAN_parse.txt',sep='\t',header=TRUE)
CAZYmarb = read.table("~/DATA/MarinobacterGenomics/2018_ProkComp/CAZY_output.txt",header=TRUE)


CAZYf  = CAZYf[!grepl("Marinobacter_", CAZYf $Genome),]


L.CAZYmarb = melt(CAZYmarb)

colnames(L.CAZYmarb) = c('CAZy_domaim', 'Genome', 'Abundance')
L.CAZYmarb = L.CAZYmarb[,c('Genome','CAZy_domaim',  'Abundance')]


CAZYf = rbind(CAZYf ,L.CAZYmarb)


unique(CAZYf$Genome)
annotatedTree$tip.label

CAZYf = CAZYf[CAZYf$Genome %in% intersect(unique(CAZYf$Genome), midRoot_bigTree$tip.label),]

CAZY = reshape(CAZYf, idvar = "Genome", timevar = "CAZy_domaim", direction = "wide")

rownames(CAZY) = CAZY[,1]
CAZY = CAZY[,2:ncol(CAZY)]
CAZY[is.na(CAZY)] <- 0
colnames(CAZY) = gsub('Abundance.','',colnames(CAZY))

tree_drop = drop.tip(midRoot_bigTree, midRoot_bigTree $tip.label[!grepl(paste(unique(CAZYf$Genome),collapse='|'), midRoot_bigTree $tip.label)])



heatmap.2(as.matrix(CAZY),trace='none',col=colorRampPalette(c('white',"#000033", "#FF3300", "#FF3300", "#FF3300"))(n = 50),margins=c(12,12))


cazyGroups = c()
for( i in c('GH','GT','PL','CE','CBM')){
	cazyGroups = cbind(cazyGroups , rowSums(CAZY[grepl(i,colnames(CAZY)),]))
}
colnames(cazyGroups) = c('GH','GT','PL','CE','CBM')

heatmap.2(as.matrix(t(cazyGroups)),trace="none", col=colorRampPalette(c('white',"#000033", "#FF3300"))(n = 50),srtCol=45,cexRow=1.5,cexCol=0.6,main='',margins=c(10,20))








genera_tip_labels = gsub("\\_.*","", tree_drop $tip.label)
genera_tip_labels = data.frame(genera_tip_labels)
colnames(genera_tip_labels)=c('Genus_name')



TaxTable$genusX1 = TaxTable$Genus_name
TaxTable$genusX2 = TaxTable$Genus_name

#Make a corrected tax thingy, as we suspect that Marinobacter as well as Marinobacterium are Oceanospirillales
TaxTableCorrected = TaxTable
levels(TaxTableCorrected $Family)=c(levels(TaxTableCorrected $Family),"Marinobacteraceae",'-')
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacter",c("Order")]= "Oceanospirillales"
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacter",c("Family")]= "Marinobacteraceae"
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacterium",c("Order")]= "Oceanospirillales"
TaxTableCorrected[TaxTableCorrected $Genus=="Marinobacterium",c("Family")]= "-"


taxa.ordered = right_join(x = TaxTableCorrected, y = genera_tip_labels, by = "Genus_name")

levels(taxa.ordered$Phylum)=c(levels(taxa.ordered$Phylum),'unknown')
levels(taxa.ordered$Class)=c(levels(taxa.ordered$Class),'unknown')
levels(taxa.ordered$Order)=c(levels(taxa.ordered$Order),'unknown')
levels(taxa.ordered$Family)=c(levels(taxa.ordered$Family),'unknown')
levels(taxa.ordered$Phylum)=c(levels(taxa.ordered$Phylum),'unknown')
levels(taxa.ordered$genusX1)=c(levels(taxa.ordered$genusX1),'unknown')
levels(taxa.ordered$genusX2)=c(levels(taxa.ordered$genusX2),'unknown')

taxa.ordered[is.na(taxa.ordered)] <- "unknown" 

rownames(taxa.ordered)= tree_drop $tip.label
taxa.ordered$strain = tree_drop $tip.label

l <- list(); 
for(Order in unique(taxa.ordered$Order)){
l[[Order]] = as.character(taxa.ordered[taxa.ordered $Order==Order,'strain'])
}

annotatedTree <- groupOTU(tree_drop, l)


p6= ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,) +geom_treescale()+theme_tree() +  scale_color_manual(values=rev(colorRampPalette(brewer.pal(8, "Set1"))(18))) + theme(legend.position="right")

t12.1 <-  gheatmap(p6, data.frame(cazyGroups), width=0.15, low="white", high="darkslategray", colnames = TRUE)
t12.1 = open_tree(t12.1, angle=180)+theme(legend.position='bottom')

t12.2 <-  gheatmap(p6, data.frame(rowSums(CAZY)), width=0.15, low="white", high="darkslategray", colnames = TRUE)



data.frame(rowSums(CAZY))





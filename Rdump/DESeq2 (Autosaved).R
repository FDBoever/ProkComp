###########################################################################################################
# 	References:
#		https://github.com/MadsAlbertsen/ampvis/blob/master/R/amp_test_species.R
#		http://userweb.eng.gla.ac.uk/umer.ijaz/
###########################################################################################################


# IF YOU NEED TO DETACH ALL PACKAGES BEFOREHAND
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

# Load the required packages
library(phyloseq)
library(randomForest)
library(dplyr)
library(ggplot2)
#library(DESeq2)

##############################################
#	COG CATEGORY LIST
##############################################

COG_meaning = data.frame(rbind(c('A','RNA processing and modification'),c('B','Chromatin Structure and dynamics'),c('C','Energy production and conversion'),c('D','Cell cycle control and mitosis'),c('E','Amino Acid metabolis and transport'),c('F','Nucleotide metabolism and transport'),c('G','Carbohydrate metabolism and transport'),c('H','Coenzyme metabolis'),c('I','Lipid metabolism'),c('J','Tranlsation'),c('K','Transcription'),c('L','Replication and repair'),c('M','Cell wall/membrane/envelop biogenesis'),c('N','Cell motility'),c('O','Post-translational modification, protein turnover, chaperone functions'),c('P','Inorganic ion transport and metabolism'),c('Q','Secondary Structure'),c('T','Signal Transduction'),c('U','Intracellular trafficing and secretion'),c('Y','Nuclear structure'),c('Z','Cytoskeleton'),c('R','General Functional Prediction only'),c('S','Function Unknown'),c('V','Defence Mechanisms'),c('W','Extracellular structures'),c('Y','Nuclear structure'),c('X','Mobilome: prophages, transposons')))
names(COG_meaning)=c("Var1","long")

##############################################


#Load in the pan-matrix (could be any nxp matrix) derived from ProkComp
panmatrix = AnvioData
meta_data =  CheckMANVIO2
read.table()

#write.table(panmatrix,"panmatrix.txt",row.names=T,col.names=T,quote=F)
#write.table(meta_data,"meta_data.txt",row.names=T,col.names=T,quote=F)
#write.table(ANVIO_cat2,"ANVIO_cat2.txt",row.names=T,col.names=T,quote=T)

panmatrix = read.table("panmatrix.txt")
meta_data = read.table("meta_data.txt")
ANVIO_cat2  = read.table("ANVIO_cat2.txt")



#trim out the unwanted (low-quality) genomes
#rownames(meta_data)=meta_data$Bin_Id
#meta_data = meta_data[,3:length(colnames(meta_data))]

#assign a new variable that contain all the grouping information
grouping_factors = ANVIO_cat2
grouping_factors = grouping_factors[grouping_factors$Bin_Id %in% rownames(meta_data),]
rownames(grouping_factors)= grouping_factors $Bin_Id
grouping_factors = grouping_factors[,2:length(colnames(grouping_factors))]

#convert the table to DESeqDataSet object
geneCount = round(as(panmatrix, "matrix"), digits = 0)

geneCountPA = geneCount 
geneCountPA[geneCountPA>1]=1


######################################################################################
#	INITIAL FILTERING STEP
#	All downstream analysis are preferably run on non-core, and non-singleton genes
#	For this purpose we select only the accessory genes to dig for patterns, and increase the statistical power
######################################################################################

CoreGenes = geneCountPA[rowSums(geneCountPA)== length(colnames(geneCountPA)), ]
none_CoreGenes = geneCountPA[!rowSums(geneCountPA)== length(colnames(geneCountPA)), ]
none_core_none_signleton = none_CoreGenes[!rowSums(none_CoreGenes)==1, ]

#geneCount only containing non-core and non-singleton genes
geneCount = geneCount[rownames(none_core_none_signleton),]

######################################################################################
#
#	Kruskal-Wallis Test with FDR
#
#	!!! calcultae which gene have differences in relative abundance (per genome) among grouping levels using non-parametric testing for >2 groups
#
#		* Using log-relative transformation or relative transformation
#				- log relative:	log((geneCount)/(rowSums(geneCount)+dim(geneCount)[2]))
#				- relative:		geneCount/rowSums(geneCount)
#		* for each gene (or) we will compute a Kruskal-Wallis test of variance
#		* p-values will be recorded and corrected using both Bonferonni and FDR corrections
#
#	Reference: 
#		http://www.bigre.ulb.ac.be/courses/statistics_bioinformatics/practicals/microarrays_berry_2010/berry_feature_selection.html
#
#######################################################################################


groups<-as.factor(grouping_factors$group2)

#### TRANSFORMATION

#using log-relative transformation
geneCount <- geneCount +1 
data<-log((geneCount)/(rowSums(geneCount)+dim(geneCount)[2]))

#using relative transformation
#data<-geneCount/rowSums(geneCount)

data<-as.data.frame(data)
data=t(data)

kruskal.wallis.alpha=0.001
kruskal.wallis.table <- data.frame()

#----------------------------------------------------------------------------------
# Calculate Kruskal-Wallis tests for each of the genes
#----------------------------------------------------------------------------------
for (i in 1:dim(data)[2]) {
 	ks.test <- kruskal.test(data[,i], g=groups)
 	kruskal.wallis.table <- rbind(kruskal.wallis.table,data.frame(id= colnames(data)[i], p.value=ks.test$p.value))
 # Report number of values tested
 cat(paste("Kruskal-Wallis test for ", colnames(data)[i]," ", i, "/", 
           dim(data)[2], "; p-value=", ks.test$p.value,"\n", sep=""))
}

#----------------------------------------------------------------------------------
#Add the corrections in additional collumns
#----------------------------------------------------------------------------------

kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, 
                                   size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor

#Alternative False Discovery Rate correction
kruskal.wallis.table$pvalsFDR = kruskal.wallis.table$p.value*(length(kruskal.wallis.table$p.value)/rank(kruskal.wallis.table$p.value,ties.method="average"))
#Bonferroni correction
kruskal.wallis.table$pvalsBon = kruskal.wallis.table$p.value*length(kruskal.wallis.table$p.value)
                                            

#Sort the dataframe according to p.value
kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value,
                                                  decreasing=FALSE), ]


#----------------------------------------------------------------------------------
#Select those genes that have FDR corrected p-values <0.001
#----------------------------------------------------------------------------------

kruskal.wallis.alpha=0.001
last.significant.element <- max(which(kruskal.wallis.table$pvalsFDR <= kruskal.wallis.alpha))
selected <- 1:last.significant.element
diff.cat.factor <- kruskal.wallis.table$id[selected]
diff.cat <- as.vector(diff.cat.factor)
Kruskal.sig = diff.cat

#Select those genes that have FDR corrected p-values <0.01
kruskal.wallis.alpha=0.01
last.significant.element <- max(which(kruskal.wallis.table$pvalsFDR <= kruskal.wallis.alpha))
selected <- 1:last.significant.element
diff.cat.factor <- kruskal.wallis.table$id[selected]
diff.cat <- as.vector(diff.cat.factor)
Kruskal.sig.0.01 = diff.cat

#Select those genes that have FDR corrected p-values <0.05
kruskal.wallis.alpha=0.05
last.significant.element <- max(which(kruskal.wallis.table$pvalsFDR <= kruskal.wallis.alpha))
selected <- 1:last.significant.element
diff.cat.factor <- kruskal.wallis.table$id[selected]
diff.cat <- as.vector(diff.cat.factor)
Kruskal.sig.0.05 = diff.cat

#Select those genes that have FDR corrected p-values <0.05
kruskal.wallis.alpha=0.000001
last.significant.element <- max(which(kruskal.wallis.table$pvalsFDR <= kruskal.wallis.alpha))
selected <- 1:last.significant.element
diff.cat.factor <- kruskal.wallis.table$id[selected]
diff.cat <- as.vector(diff.cat.factor)
Kruskal.sig.0.000001 = diff.cat


kwp <- ggplot(kruskal.wallis.table, aes(x= reorder(id,pvalsFDR), weight = pvalsFDR))
kwp <- kwp + geom_bar(color ='blue4') + ggtitle("KW-fdr")
kwp <- kwp + xlab("Protein Clusters") + ylab("p-value")
kwp <- kwp +theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
kwp <- kwp + geom_hline(aes(yintercept=0.05), colour="red", linetype="dashed")
kwp <- kwp + scale_y_continuous(expand = c(0, 0))


venn(list('0.05' = Kruskal.sig.0.05, '0.01' = Kruskal.sig.0.01,'0.001'= Kruskal.sig))


#----------------------------------------------------------------------------------
#Now we plot taxa significantly different between the categories
#----------------------------------------------------------------------------------

df<-NULL
for(i in diff.cat){
 tmp<-data.frame(data[,i],groups,rep(paste(paste(i,gsub(".*;","",gsub(";+$","",paste(collapse=";")))),"\n q = ",sprintf("%.5g",kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"]),sep=""),dim(data)[1]))
 if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","group","gene")

p<-ggplot(df[1:(57*18),],aes(group,Value,colour=group))+ylab("Log-relative normalised")
p<-p+geom_boxplot()+geom_jitter()+theme_bw()+
 facet_wrap( ~ gene )
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text.x = element_text(size = 6, colour = "black"))
p + scale_color_manual(values=c("#39811D","#86DA81","#A6D8D4","#F2A968","#F2EC70","#E38FDD","#898989","#76AECF","#B34D22"))


#----------------------------------------------------------------------------------
#Now we plot taxa significantly different between the categories
#----------------------------------------------------------------------------------
sDF4 = DF4[DF4$genes %in% Kruskal.sig.0.05,]
Kruskal.sig.0.05

COG_list = c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
	sum=0
	mat = colSums(t(table(droplevels(sDF4[grep('', sDF4$COG_CATEGORY),c('COG_CATEGORY_ACC','COG_FUNCTION')]))))

	CpcogFull = data.frame(sum(mat))
	for(i in COG_list){
		mat = colSums(t(table(droplevels(sDF4[grep(i, sDF4$COG_CATEGORY),c('COG_CATEGORY_ACC','COG_FUNCTION')]))))
		CpcogFull = cbind(CpcogFull, sum(mat))
	}
names(CpcogFull)= c('',COG_list)
CpcogFull = t(CpcogFull)
CpcogFull = cbind(CpcogFull,"COG"=c('',COG_list))
colnames(CpcogFull)=c('count','COG')
CpcogFull = data.frame(CpcogFull)
CpcogFull$count = as.numeric(as.character(CpcogFull$count ))


CpcogFull$COG2 = paste(CpcogFull$COG, COG_meaning[COG_meaning$Var1%in% CpcogFull$COG,'long'],sep=' - ')


g <- ggplot(CpcogFull, aes(x=reorder(COG2,count),y=count))
g + geom_bar(stat='identity')  + scale_y_continuous(expand = c(0, 0)) +xlab('COG category')+ coord_flip()



######################################################################################
#
#	RANDOM FOREST 
#	The Random Forests (RF) [Breiman 2001] algorithm is an increasingly popular machine learning algorithm within statistical genetics.
#	RF is among the moost accurate classification methods to date and is suitable to analyse high dimentional nxp matrices 	
#	using the randomForest() function in the 'randomForest() package
#	
#
#	Reference: 
#		doi:  10.2202/1544-6115.1691
#
#######################################################################################


set.seed(42)

#prepare the dataset suitable for randomForest()
RFdf = data.frame(t(geneCount),'group'=as.character(grouping_factors$group2))

RF <- randomForest(group ~ ., data= RFdf)
#RF <- randomForest(group ~ ., data= RFdf, importance=T, proximity=T,ntree=1500,keep.forest=F)

#TRAIN THE RF (look for code in RandomForest_Brieuc_modification.R)

RF <- randomForest(group ~ ., data= RFdf,importance=TRUE ,proximity=TRUE, ntree= 25000, mtry= 1715 ,strata= group)


# Extracts variable importance (Mean Decrease in Gini Index)
# Sorts by variable importance and relevels factors to match ordering
var_importance <- data_frame(variable=setdiff(colnames(RFdf), "group"),
                             importance=as.vector(importance(RF,type=2)))
var_importance <- arrange(var_importance, desc(importance))
var_importance$variable <- factor(var_importance$variable, levels=var_importance$variable)


p <- ggplot(var_importance, aes(x=variable, weight=importance))
p <- p + geom_bar(color="darkgreen") + ggtitle("Variable Importance from Random Forest Fit")
p <- p + xlab("Protein Clusters") + ylab("Variable Importance (Mean Decrease in Gini Index)")
p <- p +theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p <- p + geom_vline(xintercept = 1000,colour="red", linetype="dashed") + scale_y_continuous(expand = c(0, 0))

multiplot(kwp,p,cols=1)



RFselect =  var_importance[1:10,]$variable
RFselect =  var_importance[1:500,]$variable
RFselect =  as.character(var_importance[1:1000,]$variable)
#RFselect =  var_importance$variable



################################################################################################################################
#	
#	SCOARY!!!!!
#	
################################################################################################################################

##############################################
#	SCOARY INPUT GENERATOR
##############################################

library(ape)
library(phytools)

Map<-grouping_factors
Tab<-AnvioData

df_trait=cbind()
for(group in unique(Map $group2)){
	trait<-as.numeric(ifelse(Map$group2==group,1,0))
	df_trait= cbind(df_trait,trait)
}
rownames(df_trait)<-rownames(Map)
colnames(df_trait)<-unique(Map$group2)

#remove low quality genomes
df_trait = df_trait[CheckMANVIO2$Bin_Id,]

df_Tab<-data.frame(genes=rownames(Tab))
#DEDUPLICATED2 =  DEDUPLICATED[,c("protein_cluster_id",'COG_CATEGORY_ACC','COG_CATEGORY','COG_FUNCTION_ACC')]
#names(DEDUPLICATED2)= c("genes","COG_CATEGORY_ACC","COG_CATEGORY","COG_FUNCTION_ACC")
DEDUPLICATED3 =  DEDUPLICATED[,c("protein_cluster_id",'COG_CATEGORY_ACC','COG_CATEGORY','COG_FUNCTION_ACC','COG_FUNCTION')]
names(DEDUPLICATED3)= c("genes","COG_CATEGORY_ACC","COG_CATEGORY","COG_FUNCTION_ACC",'COG_FUNCTION')
DF4 = right_join(DEDUPLICATED3, df_Tab,by=('genes'))
DF4$COG_FUNCTION = sub('\\\t.*', '', DF4$COG_FUNCTION)
DF4$COG_FUNCTION = chartr("+" , "_", DF4$COG_FUNCTION)
DF4$COG_FUNCTION = chartr("," , "_", DF4$COG_FUNCTION)
DF4$COG_FUNCTION = chartr(" " , "_", DF4$COG_FUNCTION)

write.table(DF4,"DF4.txt",row.names=T,col.names=T,quote=F)

#Abundance_Table
df_Tab<-cbind(DF4,as.data.frame.matrix(Tab))

#Presence_Absence Table
Tab1 = as.data.frame.matrix(Tab)
Tab1[Tab1>0] <- 1
df_Tab_pres_abs<-cbind(DF4,Tab1)

write.table(df_trait,"traits_scoary.csv",row.names=T,col.names=NA,sep=",",append=F,quote=F)
write.table(df_Tab,"matrix_scoary_abundance.csv",row.names=F,col.names=T,sep=",",append=F,quote=F)
write.table(df_Tab_pres_abs,"matrix_scoary_Presence_Absence.csv",row.names=F,col.names=T,sep=",",append=F,quote=F)

write.table(df_Tab[rownames(geneCount),],"matrix_scoary_abundance.csv",row.names=F,col.names=T,sep=",",append=F,quote=F)

#write.table(df_Tab[,colnames(geneCount)],"matrix_scoary_abundance.tsv",row.names=F,col.names=T,sep=",",append=F,quote=F)


write.tree(drop.tip(tree2, tree2 $tip.label[-match(CheckMANVIO2$Bin_Id, tree2 $tip.label)]), file = "~/pruned_tree.newick")
write.tree(drop.tip(tree2, tree2 $tip.label[-match(colnames(geneCount), tree2 $tip.label)]), file = "~/pruned_tree")


#	RUN SCOARY IN THE TERMINAL
# /Users/sa01fd/Genomics/miniconda2/bin/scoary -t /Users/sa01fd/traits_scoary.csv -g /Users/sa01fd/matrix_scoary_Presence_Absence.csv -n /Users/sa01fd/DATA/MarinobacterGenomics/2018_ProkComp/trees/SCO_kde.fas.treefile  -e 1000 -p 1.0 --threads 8 -s 6 -o ./FINAL_scoary

#scoary -t /Users/sa01fd/traits_scoary.tsv -g /Users/sa01fd/matrix_scoary_abundance.tsv -n /Users/sa01fd/DATA/MarinobacterGenomics/2018_ProkComp/trees/SCO_kde.fas.treefile  -e 1000 -p 1.0 --threads 8 -s 6 -o ./FINAL_scoary




##############################################
#	SCOARY OUTPUT ANALYSER
##############################################
~/Genomics/multigeneblast_1.1.14_macosx_commandline/Accessory_Scoary/algicola_30_04_2018_0955.results.csv
df_CpcogFull = cbind()
Scoary.sig = c()
cladetracker = c()
plots = list()
COG_distribution_all = data.frame(cbind(Var1=as.character(unique(ANVIO$COG_CATEGORY_ACC))))
for(group in unique(Map $group2)){
	scoary <- read.table(paste('/Users/sa01fd/FINAL_scoary2/',group,'_22_06_2018_1836.results.csv',sep=''),header=TRUE,sep=",") #protein clusters
	
	Scp <- ggplot(scoary, aes(x= reorder(genes, Benjamini_H_p), weight = Benjamini_H_p))
	Scp <- Scp + geom_bar(color ='darkred') + ggtitle("")
	Scp <- Scp + xlab("Protein Clusters") + ylab("p-value")
	Scp <- Scp +theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
	
	plots[[group]]=Scp
	
	scoary = right_join(DEDUPLICATED3, scoary,by=('genes'))
	scoary = scoary[order(scoary$Benjamini_H_p),]
	scoary = scoary[scoary$Benjamini_H_p <= 0.05,]
	scoary = scoary[scoary$Number_pos_present_in > scoary$Number_pos_not_present_in,]

	COG_dist = data.frame(table(scoary $COG_CATEGORY_ACC.x))
	COG_dist  = COG_dist[COG_dist$Freq>0,]
	COG_dist = right_join(COG_meaning,COG_dist,by=('Var1'))
	COG_dist$combined = paste(COG_dist$Var1, COG_dist$long)
	p<-ggplot(COG_dist, aes(reorder(combined,Freq), Freq)) + geom_bar(stat = "identity") +theme(axis.text.x=element_text(angle=90,hjust=1))+coord_flip()+ggtitle(group)
	p
	COG_distribution_all = right_join(COG_dist[,c('Var1','Freq')], COG_distribution_all,by=('Var1'))


	COG_list = c("D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
	sum=0
	mat = colSums(t(table(droplevels(scoary[grep("C", scoary$COG_CATEGORY.y),c('COG_CATEGORY_ACC.x','COG_FUNCTION')]))))
	CpcogFull = data.frame(sum(mat))
	for(i in COG_list){
		mat = colSums(t(table(droplevels(scoary[grep(i, scoary$COG_CATEGORY.y),c('COG_CATEGORY_ACC.x','COG_FUNCTION')]))))
		CpcogFull = cbind(CpcogFull, sum(mat))
	}
	df_CpcogFull = rbind(df_CpcogFull, CpcogFull)
	Scoary.sig = c(Scoary.sig,as.character(scoary$genes))
	cladetracker = c(cladetracker, rep(group,length(scoary$genes)))	
}

names(df_CpcogFull)= c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
rownames(df_CpcogFull)=as.character(unique(ANVIO_cat$group2))
df_CpcogFull = data.frame(t(df_CpcogFull))
df_CpcogFull$Var1 = c( "C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )

df_CpcogFull = right_join(COG_meaning, df_CpcogFull,by=('Var1'))
df_CpcogFull2 = df_CpcogFull[,3:ncol(df_CpcogFull)]
rownames(df_CpcogFull2) = paste(df_CpcogFull$Var1 , df_CpcogFull$long, sep=' - ')

heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=viridis(50),margins=c(12,12), srtCol =45,ColSideColors=unique(cladecolors))
heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=colorRampPalette(c("#000033", "#FF3300"))(n = 50),margins=c(12,12), srtCol =45,ColSideColors=unique(cladecolors))
heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=colorRampPalette(c("#F7FBFF", "#000033"))(n = 50),margins=c(12,12), srtCol =45,ColSideColors=unique(cladecolors))


ScoaryEnrichedGenes = data.frame('gene' = Scoary.sig, 'clade' = cladetracker,'E_D' = rep('Enriched',length(Scoary.sig)))



######### ---- HEATMAP ---- Enriched COGs #########


ANVIOfiltered = ANVIO[ANVIO$genome_name %in% CheckMANVIO2$Bin_Id,]
GroupB <-as.character(scoary$genes)

mat = t(table(droplevels(subset(ANVIOfiltered, protein_cluster_id  %in% ScoaryEnrichedGenes$gene)[,c('COG_FUNCTION','genome_name')])))
pruned_tree = drop.tip(tree2, tree2 $tip.label[-match(CheckMANVIO2$Bin_Id, tree2 $tip.label)])
pruned_tree $edge.length[which(pruned_tree $edge.length == 0)] <- 0.00001
pruned_tree_um <- chronopl(pruned_tree,lambda = 0.1,tol = 0)
pruned_tree_d <- as.dendrogram(as.hclust.phylo(pruned_tree_um))

clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

heatmap.2(as.matrix(t(combined_ordered_matrix)),trace="none"
,col=colorRampPalette(c("white","black","black","black","black","black"))(n = 150),srtCol=45,cexRow=0.6,cexCol=0.6,main='',margins=c(10,20),ColSideColors= cladecolors[new_order],Colv= pruned_tree_d)


#######################
#	ANALYLISIS DEPLETED
#######################

df_CpcogFull = cbind()
Scoary.sig = c()
cladetracker = c()
COG_distribution_all = data.frame(cbind(Var1=as.character(unique(ANVIO$COG_CATEGORY_ACC))))
for(group in unique(Map $group2)){
	scoary <- read.table(paste('/Users/sa01fd/DATA/MarinobacterGenomics/2018_ProkComp/Scoary/',group,'_24_04_2018_2026.results.csv',sep=''),header=TRUE,sep=",") #protein clusters
	scoary = right_join(DEDUPLICATED3, scoary,by=('genes'))
	scoary = scoary[order(scoary$Benjamini_H_p),]
	scoary = scoary[scoary$Benjamini_H_p <= 0.01,]
	
	#Crucual line to change
	scoary = scoary[scoary$Number_pos_present_in < scoary$Number_pos_not_present_in,]

	COG_dist = data.frame(table(scoary $COG_CATEGORY_ACC.x))
	COG_dist  = COG_dist[COG_dist$Freq>0,]
	COG_dist = right_join(COG_meaning,COG_dist,by=('Var1'))
	COG_dist$combined = paste(COG_dist$Var1, COG_dist$long)
	p<-ggplot(COG_dist, aes(reorder(combined,Freq), Freq)) + geom_bar(stat = "identity") +theme(axis.text.x=element_text(angle=90,hjust=1))+coord_flip()+ggtitle(group)
	p
	COG_distribution_all = right_join(COG_dist[,c('Var1','Freq')], COG_distribution_all,by=('Var1'))


	COG_list = c("D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
	sum=0
	mat = colSums(t(table(droplevels(scoary[grep("C", scoary$COG_CATEGORY.y),c('COG_CATEGORY_ACC.x','COG_FUNCTION')]))))
	CpcogFull = data.frame(sum(mat))
	for(i in COG_list){
		mat = colSums(t(table(droplevels(scoary[grep(i, scoary$COG_CATEGORY.y),c('COG_CATEGORY_ACC.x','COG_FUNCTION')]))))
		CpcogFull = cbind(CpcogFull, sum(mat))
	}
	df_CpcogFull = rbind(df_CpcogFull, CpcogFull)
	Scoary.sig = c(Scoary.sig,as.character(scoary$genes))
	cladetracker = c(cladetracker, rep(group,length(scoary$genes)))

}

names(df_CpcogFull)= c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
rownames(df_CpcogFull)=as.character(unique(ANVIO_cat$group2))
df_CpcogFull = data.frame(t(df_CpcogFull))
df_CpcogFull$Var1 = c( "C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )

df_CpcogFull = right_join(COG_meaning, df_CpcogFull,by=('Var1'))
df_CpcogFull2 = df_CpcogFull[,3:ncol(df_CpcogFull)]
rownames(df_CpcogFull2) = paste(df_CpcogFull$Var1 , df_CpcogFull$long, sep=' - ')

heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=viridis(50),margins=c(12,12), srtCol =45,ColSideColors=unique(cladecolors))
heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=colorRampPalette(c("#FFF5F0", "#67000D"))(n = 50),margins=c(12,12), srtCol =45,cexRow = 0.8, cexCol = 0.8)
heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=colorRampPalette(c("FFF5F0", "#67000D"))(n = 50),margins=c(12,12), srtCol =45,ColSideColors=unique(cladecolors))

ScoaryDepletedGenes = data.frame('gene' = Scoary.sig, 'clade' = cladetracker,'E_D' = rep('Depleted',length(Scoary.sig)))

######### ---- HEATMAP ---- Enriched COGs #########

cladecolorsHeat = cladecolors
names(cladecolorsHeat)=CheckMANVIO2$Bin_Id

mat = t(table(droplevels(subset(ANVIOfiltered, protein_cluster_id  %in% ScoaryDepletedGenes $gene)[,c('COG_FUNCTION','genome_name')])))
pruned_tree = drop.tip(tree2, tree2 $tip.label[-match(CheckMANVIO2$Bin_Id, tree2 $tip.label)])
pruned_tree $edge.length[which(pruned_tree $edge.length == 0)] <- 0.00001
pruned_tree_um <- chronopl(pruned_tree,lambda = 0.1,tol = 0)
pruned_tree_d <- as.dendrogram(as.hclust.phylo(pruned_tree_um))

clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

heatmap.2(as.matrix(t(combined_ordered_matrix)),trace="none"
,col=colorRampPalette(c("white","black","black","black","black","black"))(n = 150),srtCol=45,cexRow=0.6,cexCol=0.6,main='',margins=c(10,20),ColSideColors= as.character(cladecolorsHeat[rownames(combined_ordered_matrix)]),Colv= pruned_tree_d)

######### ---- HEATMAP ---- BOTH ENRICHED AND DEPLETED COGs #########

cladecolorsHeat = cladecolors
names(cladecolorsHeat)=CheckMANVIO2$Bin_Id
meta_data$group2 = as.character(ANVIO_cat[ANVIO_cat$name %in% CheckMANVIO2$Bin_Id,]$group2)




mat = t(table(droplevels(subset(ANVIOfiltered, protein_cluster_id  %in% rbind(ScoaryEnrichedGenes, ScoaryDepletedGenes)$gene)[,c('COG_FUNCTION','genome_name')])))
pruned_tree = drop.tip(tree2, tree2 $tip.label[-match(CheckMANVIO2$Bin_Id, tree2 $tip.label)])
pruned_tree $edge.length[which(pruned_tree $edge.length == 0)] <- 0.00001
pruned_tree_um <- chronopl(pruned_tree,lambda = 0.1,tol = 0)
pruned_tree_d <- as.dendrogram(as.hclust.phylo(pruned_tree_um))

clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

heatmap.2(as.matrix(t(combined_ordered_matrix)),trace="none"
,col=colorRampPalette(c("white","black","black","black","black","black"))(n = 150),srtCol=45,cexRow=0.6,cexCol=0.6,main='',margins=c(10,20),ColSideColors= as.character(cladecolorsHeat[rownames(combined_ordered_matrix)]),Colv= pruned_tree_d)

#### -----
library(ggridges)

ScoaryEnrichedDepleted = rbind(ScoaryEnrichedGenes, ScoaryDepletedGenes)

ScoaryEnrichedDepleted

 meta_data$group2 = as.character(ANVIO_cat[ANVIO_cat$name %in% CheckMANVIO2$Bin_Id,]$group2)
 
#ggplot(meta_data, aes(x = GC, y = group2)) +
#  geom_density_ridges(aes(fill = group2)) +
#  scale_fill_manual(values = unique(cladecolors))

scp = ggplot(ScoaryEnrichedDepleted, aes(E_D, ..count..)) + geom_bar(aes(fill = clade), position = "dodge")+scale_fill_manual(values = c("#39811D","#86DA81","#A6D8D4","#F2A968","#F2EC70", '#888888',"#76AECF","#B34D22"))+ scale_y_continuous(expand = c(0, 0))

multiplot(kwp,p, scp)


#ggplot(ScoaryEnrichedDepleted, aes(x = E_D, y = ..count..)) +
#  geom_density(aes(fill = clade)) +
#  scale_fill_manual(values = unique(cladecolors))




##################################################################

##################################################################


venn(list('KruskalFDR' = Kruskal.sig, 'Random Forest' = RFselect))
venn(list('KW+FDR<0.05' = Kruskal.sig.0.05, 'KW+FDR<0.01' = Kruskal.sig.0.01,'KW+FDR<0.001'= Kruskal.sig,'RandomForest[2000]'=RFselect))
venn(list('KW+FDR<0.05' = Kruskal.sig.0.05, 'KW+FDR<0.01' = Kruskal.sig.0.01,'KW+FDR<0.001'= Kruskal.sig,'RandomForest[2000]'=RFselect,'Scoary'= Scoary.sig))

Scoary.sig = intersect(Scoary.sig ,rownames(geneCount))

intersect(Scoary.sig , RFselect, Kruskal.sig.0.05)

mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in% RFselect)[,c('COG_FUNCTION','genome_name')])))


mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in% RFselect)[,c('protein_cluster_id','genome_name')])))
mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in% intersect(Kruskal.sig, RFselect))[,c('protein_cluster_id','genome_name')])))
mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in% intersect(kruskal.wallis.table[kruskal.wallis.table$q.value < 0.05 ,c('id')], RFselect))[,c('protein_cluster_id','genome_name')])))

mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in% intersect(kruskal.wallis.table[kruskal.wallis.table$q.value < 0.01 ,c('id')], RFselect))[,c('protein_cluster_id','genome_name')])))

mat = mat[CheckMANVIO2$Bin_Id,]

mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in% Scoary.sig)[,c('protein_cluster_id','genome_name')])))
mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in%  intersect(Scoary.sig,rownames(geneCount)))[,c('protein_cluster_id','genome_name')])))
mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in%  intersect(RFselect ,intersect(Scoary.sig,rownames(geneCount))))[,c('protein_cluster_id','genome_name')])))


Kruskal.sig.0.000001

mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in% Kruskal.sig.0.000001)[,c('COG_FUNCTION','genome_name')])))

mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in% Kruskal.sig)[,c('COG_FUNCTION','genome_name')])))
pruned_tree = drop.tip(tree2, tree2 $tip.label[-match(CheckMANVIO2$Bin_Id, tree2 $tip.label)])
pruned_tree $edge.length[which(pruned_tree $edge.length == 0)] <- 0.00001
pruned_tree_um <- chronopl(pruned_tree,lambda = 0.1,tol = 0)
pruned_tree_d <- as.dendrogram(as.hclust.phylo(pruned_tree_um))

clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

heatmap.2(as.matrix(t(combined_ordered_matrix)),trace="none"
,col=colorRampPalette(c("white","black"))(n = 150),srtCol=45,cexRow=0.6,cexCol=0.6,main='',margins=c(10,20),ColSideColors= cladecolors[new_order],Colv= pruned_tree_d)
heatmap.2(as.matrix(t(combined_ordered_matrix)),trace="none"
,col=colorRampPalette(c("white","black",'black','#000033','#FF3300'))(n = 150),srtCol=45,cexRow=0.6,cexCol=0.6,main='',margins=c(10,20),ColSideColors= cladecolors[new_order],Colv= pruned_tree_d)

heatmap.2(as.matrix(log(combined_ordered_matrix+1)),trace="none"
,col=colorRampPalette(c("#30357B",'white','#DD8920','#DD8920','#DD8920'))(n = 150),srtCol=45,cexRow=0.6,cexCol=0.6,main='',margins=c(10,20),RowSideColors= cladecolors[new_order],Rowv= pruned_tree_d)


heatmap.2(as.matrix(t(combined_ordered_matrix)),trace="none"
,col = inferno(75),srtCol=45,cexRow=0.6,cexCol=0.6,main='',margins=c(10,20),ColSideColors= cladecolors[new_order],Colv= pruned_tree_d)

heatmap.2(as.matrix(t(combined_ordered_matrix)),trace="none"
,col = inferno(75),srtCol=45,cexRow=0.6,cexCol=0.6,main='',margins=c(10,20),ColSideColors= cladecolors[new_order],Colv= pruned_tree_d)

"white","black",'black','#000033',#FF3300"


#########
intersect = intersect( as.character(RFselect), as.character(ScoaryEnrichedDepleted$gene))
intersect = intersect(Kruskal.sig.0.05, intersect)


mat = t(table(droplevels(subset(ANVIOfiltered, protein_cluster_id  %in% intersect)[,c('COG_FUNCTION','genome_name')])))
pruned_tree = drop.tip(tree2, tree2 $tip.label[-match(CheckMANVIO2$Bin_Id, tree2 $tip.label)])
pruned_tree $edge.length[which(pruned_tree $edge.length == 0)] <- 0.00001
pruned_tree_um <- chronopl(pruned_tree,lambda = 0.1,tol = 0)
pruned_tree_d <- as.dendrogram(as.hclust.phylo(pruned_tree_um))

clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

heatmap.2(as.matrix(t(combined_ordered_matrix)),trace="none"
,col=colorRampPalette(c("white","black","black","black","black","black"))(n = 150),srtCol=45,cexRow=0.6,cexCol=0.6,main='',margins=c(10,20),ColSideColors= as.character(cladecolorsHeat[rownames(combined_ordered_matrix)]),Colv= pruned_tree_d)


##########


mat = t(table(droplevels(subset(ANVIOfiltered, protein_cluster_id  %in% intersect)[,c('protein_cluster_id','genome_name')])))
pruned_tree = drop.tip(tree2, tree2 $tip.label[-match(CheckMANVIO2$Bin_Id, tree2 $tip.label)])
pruned_tree $edge.length[which(pruned_tree $edge.length == 0)] <- 0.00001
pruned_tree_um <- chronopl(pruned_tree,lambda = 0.1,tol = 0)
pruned_tree_d <- as.dendrogram(as.hclust.phylo(pruned_tree_um))

clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]


RFselect =  as.character(var_importance[1:1000,]$variable)
all.sig  = unique(c(as.character(RFselect), as.character(ScoaryEnrichedDepleted$gene), as.character(Kruskal.sig.0.05)))

DF_significance = c()
for(i in unique(ANVIOfiltered$protein_cluster_id)){
	DF_significance = rbind(DF_significance ,c(i , i %in% RFselect, i %in% ScoaryEnrichedDepleted$gene, i %in% Kruskal.sig.0.05))	
}
colnames(DF_significance) = c('gene','RF','Scoary','KruskalFDR')

DF_sig = DF_significance
colnames(DF_sig) = c('genes','RF','Scoary','KruskalFDR')
#Set Colors for RandomForest Significance
DF_significance[,2][DF_significance[,2]=='TRUE'] = 'darkgreen'
DF_significance[,2][DF_significance[,2]=='FALSE'] = 'darkseagreen1'

#Set Colors for Scoary Significance
DF_significance[,3][DF_significance[,3]=='TRUE'] = 'darkred'
DF_significance[,3][DF_significance[,3]=='FALSE'] = 'pink'

#Set Colors for KruskalFDR Significance
DF_significance[,4][DF_significance[,4]=='TRUE'] = 'blue4'
DF_significance[,4][DF_significance[,4]=='FALSE'] = 'lightskyblue1'


DF_significance = data.frame(DF_significance)

rownames(DF_significance) = DF_significance $gene
DF_significance = DF_significance[,2:length(colnames(DF_significance))]

DF_sig = right_join(data.frame(DF_sig), DEDUPLICATED3, by = 'genes')



COGselect = DF_sig[DF_sig$COG_CATEGORY =='G',]
COGselect = DF_sig[DF_sig$COG_CATEGORY =='E',]

intersect(all.sig, COGselect[,'genes'])

COGselect[COGselect$genes %in% all.sig,'genes']

mat = t(table(droplevels(subset(ANVIOfiltered, protein_cluster_id  %in%  COGselect[COGselect$genes %in% all.sig,'genes']  )[,c('COG_FUNCTION','genome_name')])))
pruned_tree = drop.tip(tree2, tree2 $tip.label[-match(CheckMANVIO2$Bin_Id, tree2 $tip.label)])
pruned_tree $edge.length[which(pruned_tree $edge.length == 0)] <- 0.00001
pruned_tree_um <- chronopl(pruned_tree,lambda = 0.1,tol = 0)
pruned_tree_d <- as.dendrogram(as.hclust.phylo(pruned_tree_um))


for(i in c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X")){
	COGselect = DF_sig[DF_sig$COG_CATEGORY == i ,]
    pdf(paste("PC_all.sig","_", i , ".pdf", sep = ""))

	heatmap.plus(as.matrix((combined_ordered_matrix[COGselect[COGselect$genes %in% all.sig,'genes'],])), col=colorRampPalette(c('white','black'))(n = 150),Colv= pruned_tree_d,RowSideColors= as.matrix(DF_significance[COGselect[COGselect$genes %in% all.sig,'genes'],]),cexRow=0.6,cexCol=0.6,labRow=COGselect[COGselect$genes %in% all.sig,'COG_FUNCTION'],scale='none',margins=c(10,20),main=paste(COG_meaning[COG_meaning$Var1 == i,1],COG_meaning[COG_meaning$Var1 == i,2],sep=' - '))
		dev.off()
		
	pdf(paste("COG_grouped_all.sig","_", i , ".pdf", sep = ""))

	mat = t(table(droplevels(subset(ANVIOfiltered, protein_cluster_id  %in%  COGselect[COGselect$genes %in% all.sig,'genes']  )[,c('COG_FUNCTION','genome_name')])))
	
	heatmap.plus(as.matrix(t(mat)), col=colorRampPalette(c('white','black'))(n = 150),cexRow=0.6,cexCol=0.6,scale='none',margins=c(10,20),main=paste(COG_meaning[COG_meaning$Var1 == i,1],COG_meaning[COG_meaning$Var1 == i,2],sep=' - '))
	
	
	dev.off()

}

heatmap.plus(as.matrix((combined_ordered_matrix[COGselect[COGselect$genes %in% all.sig,'genes'],])), col=colorRampPalette(c('white','black'))(n = 150),Colv= pruned_tree_d,RowSideColors= as.matrix(DF_significance[COGselect[COGselect$genes %in% all.sig,'genes'],]),cexRow=0.6,cexCol=0.6,labRow=COGselect[COGselect$genes %in% all.sig,'COG_FUNCTION'],scale='none',margins=c(10,20),main='E')

cladecolorsHeat2 = as.matrix(cladecolorsHeat)
cladecolorsHeat2 = data.frame(cladecolorsHeat2)

geneSelector = "TRAP"
geneSelector = "mannitol/chloroaromatic"
geneSelector = "permease"
geneSelector = "sugar"
geneSelector = "phosphotransferase "
geneSelector = "proline "
geneSelector = "serine "

mat = t(table(droplevels(ANVIOfiltered[grepl(geneSelector, ANVIOfiltered$COG_FUNCTION),][,c('COG_FUNCTION','genome_name')])))
clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

	heatmap.plus(as.matrix(t(combined_ordered_matrix)), col=colorRampPalette(c('white','black'))(n = 150),cexRow=0.6,cexCol=0.6,scale='none',margins=c(10,20),Colv= pruned_tree_d)
dim(mat)

mat = t(table(droplevels(ANVIOfiltered[grepl(geneSelector, ANVIOfiltered$COG_FUNCTION),][,c('protein_cluster_id','genome_name')])))
clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

	heatmap.plus(as.matrix(t(combined_ordered_matrix)), col=colorRampPalette(c('white','black'))(n = 150),cexRow=0.6,cexCol=0.6,scale='none',margins=c(10,20),Colv= pruned_tree_d)

dim(mat)


sequenceSelector ='MMSSGLLALLAFTPILLAGVLLIGMRWPARRAMPVVFLVTALIGYSAWDMTLNRILASTLQGLVITIGLLWIIFGAILLLNTLKHSGAITTIRAGFTNITPDRRIQAIIIVWLFGSFIEGASGFGTPAAIAAPLLVAVGFPAMAAVMLGMMVQSTPVSFGAVGTPIVVGVTTGLDRASITERLEALGSNWDTYLQLITSEVAITHAIVGVVMPVLMVTMMTRFFGKNKSWKEGLEVLPFALFAGIAFVVPYALTGVILGPEFPSLLGGLIGLAIVTTAAKKGFLIPKRTWDFADSKDWPSEWLGSIEMKLEDIAAKPMSGFRAWVPYVLVGVVLVLSRTVEPIKQAFTSVGVSFSNILGEAGINAGIQPLYLPGGILVMVVLATFFIHRMNLKALNSAVKESGGVLLSAGFVLLFTVPMVRILINSGVNLSDLPSMPLAMATWAADAVGGIYPLLAPTVGALGAFLAGSNTVSNMMFSQFQFGVAESLGLSTALMVAVQAVGAAAGNMVAIHNVVAASATVGLLGREGQTLRKTVWPTLYYLVMTGSIALLAAYGLGLTDPLLN'

ANVIOfiltered[grepl(sequenceSelector, ANVIOfiltered$COG_FUNCTION),]



mat = t(table(droplevels(ANVIOfiltered[grepl('COG1250', ANVIOfiltered$COG_FUNCTION_ACC),][,c('COG_FUNCTION','genome_name')])))
clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

	heatmap.plus(as.matrix(t(combined_ordered_matrix)), col=colorRampPalette(c('white','black'))(n = 150),cexRow=0.6,cexCol=0.6,scale='none',margins=c(10,20),Colv= pruned_tree_d)

dim(mat)

mat = t(PfamData[grepl('perm',rownames(PfamData)),])
clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

	heatmap.plus(as.matrix(t(combined_ordered_matrix)), col=colorRampPalette(c('white','black'))(n = 150),cexRow=0.6,cexCol=0.6,scale='none',margins=c(10,20),Colv= pruned_tree_d)



AA_permease_2 
Trp_Tyr_perm
GntP_permease
AA_permease
Xan_ur_permease
Lactate_perm


toMatch <- c("AA_permease_2", "Trp_Tyr_perm", "GntP_permease",'AA_permease','Xan_ur_permease','Lactate_perm')
toMatch <- c('3HCDH﻿', '3HCDH_N﻿', 'ECH_1','Thiolase_C','Thiolase_N','FMO','FA_desaturase','Abhydrolase')
toMatch=c("Photo_RC","APS_kinase","ATP_sulfurylase","BChl_A","RuBisCO_large","RuBisCO_small","COXG","Nitr_red_alph_N","Nitrate_red_gam","NosL","Cytochrom_C552","APS_reductase_C","NIR_SIR","DmsC","PFO_beta_C","RC_P840_PscD","Citrate_synt","MeMo_Hyd_G","MCR_alpha","MCR_alpha_N","MCR_beta","MCR_beta_N","MCR_gamma","AMO","AmoA","Archaeal_AmoA","AmoC","Monooxygenase_B","NIR_SIR_ferr","Gln_synt_C","Gln_synt_N","Glu_syn_central","Glu_synthase","Fer4_NifH","Nitro_FeMo_Co","Bac_GDH","GDH_N","NAD_GH","COX1","Cytochrom_C","PRK","COX3","Cytochrome_CBB3","CitF","CDO_I","CdhC")

toMatch=c('Cytochr','COX')
toMatch=c('Glycos_transf_2','Glycos_transf_1','SLT','Transpeptidase','Glyco_transf_5','CBM_48','CBM_48','Glyco_hydro_2_N','Esterase','Glyco_hydro_3_C','Alpha-amylase_C','Chitin_synth_2','Glyco_hydro_2_C','Chitin_synth_1','Glyco_hydro_43','Cellulose_synt','Glycos_transf_N','Amidohydro_1','Glyco_hydro_32C','3D','Bgal_small_N','Glyco_hydro_38C','Glyco_hydro_65N','Glyco_hydro_92','XET_C','SBP_bac_3','Alpha-L-AF_C','Abhydrolase_3','Rod-binding','BiPBP_C','Alpha-mann_mid','Glyco_hydro_42M','PG_binding_1','MLTD_N','Esterase_phd','Chitin_bind_1','DUF1205','Alpha-amylase_N','Gly_transf_sug','DUF847','Trehalose_PPase','Branch','COesterase','GT36_AF','MGDG_synth','Alg14','Lyase_8_N','DUF1972','Sucrose_synth','Glyco_hydro_20b','DUF1975','PilZ','X8','Glyco_transf_8N','DUF1957','DUF2029','Glyco_hydro_65C','Bac_rhamnosid_N','PG_binding_3','Glyco_hydro_53','ChitinaseA_N','PUD','Fuc4NAc_transf','Glyco_hydro_97','CBM_X','PUD','AXE1','MIR','Glyco_transf_36','Caps_synth','PUD','Chitin_synth_1N','Lyase_8_C','AMPKBI','GtrA','NodZ','Glycogen_syn','CHB_HEX','Alginate_lyase2','NAGidase','Gal_Lectin','Glyco_hydro_42C','Glyco_hydro_67M','Glyco_hydro_67C','CST-I','Raffinose_syn','Pec_lyase','CelD_N','Glyco_transf_52','Glucodextran_N','CHB_HEX_C','GFO_IDH_MocA','S6PP','DUF1083','Mannosyl_trans2','DUF563','PAE','DUF821','PMEI','Sialidase','CgtA','RGP','DUF303','Hyaluronidase_1')


mat = t(PfamData[grepl(paste(toMatch,collapse="|"),rownames(PfamData)),])
clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

	heatmap.plus(as.matrix(t(combined_ordered_matrix)), col=colorRampPalette(c('white','black'))(n = 150),cexRow=0.6,cexCol=0.6,scale='none',margins=c(10,20),Colv= pruned_tree_d)
'Glycos_transf_2 ','Glycos_transf_1 ','SLT ','Transpeptidase ','Glyco_transf_5 ','CBM_48 ','CBM_48 ','Glyco_hydro_2_N ','Esterase ','Glyco_hydro_3_C ','Alpha-amylase_C ','Chitin_synth_2 ','Glyco_hydro_2_C ','Chitin_synth_1 ','Glyco_hydro_43 ','Cellulose_synt ','Glycos_transf_N ','Amidohydro_1 ','Glyco_hydro_32C ','3D ','Bgal_small_N ','Glyco_hydro_38C ','Glyco_hydro_65N ','Glyco_hydro_92 ','XET_C ','SBP_bac_3 ','Alpha-L-AF_C ','Abhydrolase_3 ','Rod-binding ','BiPBP_C ','Alpha-mann_mid ','Glyco_hydro_42M ','PG_binding_1 ','MLTD_N ','Esterase_phd ','Chitin_bind_1 ','DUF1205 ','Alpha-amylase_N ','Gly_transf_sug ','DUF847 ','Trehalose_PPase ','Branch ','COesterase ','GT36_AF ','MGDG_synth ','Alg14 ','Lyase_8_N ','DUF1972 ','Sucrose_synth ','Glyco_hydro_20b ','DUF1975 ','PilZ ','X8 ','Glyco_transf_8N ','DUF1957 ','DUF2029 ','Glyco_hydro_65C ','Bac_rhamnosid_N ','PG_binding_3 ','Glyco_hydro_53 ','ChitinaseA_N ','PUD ','Fuc4NAc_transf ','Glyco_hydro_97 ','CBM_X ','PUD ','AXE1 ','MIR ','Glyco_transf_36 ','Caps_synth ','PUD ','Chitin_synth_1N ','Lyase_8_C ','AMPKBI ','GtrA ','NodZ ','Glycogen_syn ','CHB_HEX ','Alginate_lyase2 ','NAGidase ','Gal_Lectin ','Glyco_hydro_42C ','Glyco_hydro_67M ','Glyco_hydro_67C ','CST-I ','Raffinose_syn ','Pec_lyase ','CelD_N ','Glyco_transf_52 ','Glucodextran_N ','CHB_HEX_C ','GFO_IDH_MocA ','S6PP ','DUF1083 ','Mannosyl_trans2 ','DUF563 ','PAE ','DUF821 ','PMEI ','Sialidase ','CgtA ','RGP ','DUF303 ','Hyaluronidase_1'


require(data.table)
PfamConvert <-as.data.frame(fread("~/DATA/MarinobacterGenomics/2018_ProkComp/Pfam-A.clans.tsv"))
colnames(PfamConvert) = c('ID','Clan','ac1','acc','Function')
rownames(PfamConvert) = PfamConvert$ID

Pfam2GO <-as.data.frame(fread("~/DATA/MarinobacterGenomics/2018_ProkComp/pfam2go.txt",header=FALSE))
colnames(Pfam2GO) = c('ID','acc','GOfunc','GO')

Pfam2GO[grepl('cobalamin',Pfam2GO$GOfunc,ignore.case=TRUE),'ID']



SulfurCycle =c('PF00890','PF02910','PF12139','PF00174','PF03404','PF01087','PF01747','PF14306','PF12838','PF00581','PF07680','PF04173','PF00034','PF02872','PF13501','PF08770','PF13442','PF00581','PF13442','PF01011','PF13360','PF01077','PF03460','PF14697','PF00384','PF01568','PF04879','PF12800','PF13247','PF01292','PF00384','PF01568','PF04879','PF13247','PF03916','PF17179','PF00175','PF10418','PF13247','PF03916','PF00890','PF02910','PF12139','PF01077','PF03460','PF03460','PF00037','PF01077','PF03460','PF12831','PF02662','PF12831','PF02665','PF13187','PF02665','PF04358','PF13183','PF12838','PF00384','PF01568','PF04879','PF10518','PF13247','PF03916','PF00384','PF13247','PF01568','PF03916','PF10518','PF04879','PF00374','PF00374','PF01058','PF01058','PF00175','PF00970','PF10418','PF17179','PF00581','PF00581','PF00384','PF01568','PF04879','PF12800','PF13247','PF01292','PF00581','PF00484','PF03358','PF00296','PF09084','PF00528','PF00005','PF00175','PF00258','PF00667','PF01077','PF03460','PF01507','PF01077','PF03460','PF01077','PF03460','PF01583','PF00009','PF01507','PF00459','PF00005','PF12857','PF00528','PF00528','PF13531','PF01384','PF13247','PF03916','PF00890','PF02910','PF12139','PF00174','PF03404','PF02635','PF02635','PF00556','PF00556','PF00124','PF00124','PF03460','PF01077','PF03460','PF01077','PF03460','PF03967','PF05239','PF11511','PF07992','PF07992','PF09242','PF10518','PF05398','PF13183','PF12838','PF04077','PF04295','PF00171','PF00296','PF01613','PF00501','PF13193','PF02775','PF02776','PF02615','PF02775','PF01515','PF01515','PF02776','PF00389','PF02826','PF02615','PF04295','PF00291','PF00205','PF01070','PF00296','PF00171','PF00384','PF01568','PF13247','PF09459','PF00296','PF01613','PF00296','PF00175','PF00970','PF02332','PF02406','PF04663','PF04945','PF06099','PF00378','PF00501','PF13193','PF00557','PF01571','PF08669','PF02515','PF16867','PF00441','PF02770','PF02771','PF12806','PF00378','PF00501','PF13193','PF00441','PF02770','PF02771','PF12806','PF01515','PF01515','PF00106','PF00171','PF00205','PF02775','PF02776','PF02635','PF02635','PF03460','PF01077','PF03460','PF13183','PF04077','PF07682','PF00355','PF00111','PF00175','PF00970','PF00355','PF00848','PF13577','PF02668','PF00202','PF01515','PF01515','PF00005','PF00202','PF04069','PF00528','PF00205','PF02775','PF02776','PF00171','PF01208','PF01208','PF01208','PF01370','PF00534','PF01553','PF13439','PF00374','PF01058','PF01568','PF01747','PF02662','PF04358','PF04879','PF08770','PF09242','PF10418','PF12139','PF13247')

Oxygen =c('PF00067','PF00115','PF00116','PF00117','PF00120','PF00142','PF00148','PF00173','PF00174','PF00175','PF00194','PF00208','PF00310','PF00384','PF00441','PF00484','PF00696','PF00795','PF00970','PF00988','PF01077','PF01493','PF01568','PF01645','PF02142','PF02239','PF02335','PF02461','PF02560','PF02665','PF02771','PF02786','PF02787','PF02812','PF03060','PF03063','PF03069','PF03139','PF03404','PF03460','PF03951','PF04324','PF04879','PF04898','PF05088','PF07732','PF07992','PF11844','PF13247','PF13442','PF13447','PF13806','PF14691','PF14710','PF14711')

Nitrogen =c('PF00034','PF00067','PF00115','PF00116','PF00120','PF00142','PF00146','PF00148','PF00173','PF00174','PF00175','PF00202','PF00205','PF00208','PF00258','PF00310','PF00329','PF00346','PF00355','PF00361','PF00384','PF00420','PF00441','PF00484','PF00499','PF00507','PF00662','PF00696','PF00795','PF00890','PF00970','PF01035','PF01058','PF01077','PF01257','PF01266','PF01292','PF01493','PF01512','PF01515','PF01568','PF01645','PF02239','PF02335','PF02461','PF02560','PF02665','PF02668','PF02754','PF02771','PF02775','PF02776','PF02812','PF03060','PF03063','PF03069','PF03139','PF03264','PF03404','PF03460','PF03892','PF03951','PF04324','PF04744','PF04879','PF04896','PF04898','PF05088','PF07085','PF07732','PF07992','PF09163','PF10531','PF10588','PF10589','PF11844','PF12437','PF12680','PF12800','PF12801','PF12838','PF12942','PF13183','PF13237','PF13247','PF13435','PF13442','PF13447','PF13500','PF13510','PF13806','PF14691','PF14710','PF14711','PF16694','PF16901')
Carbon =c('PF00037','PF00067','PF00115','PF00116','PF00117','PF00120','PF00142','PF00146','PF00148','PF00173','PF00174','PF00175','PF00194','PF00208','PF00296','PF00310','PF00329','PF00346','PF00361','PF00374','PF00384','PF00420','PF00441','PF00484','PF00499','PF00507','PF00662','PF00696','PF00795','PF00871','PF00970','PF00988','PF01011','PF01058','PF01077','PF01208','PF01266','PF01493','PF01515','PF01568','PF01571','PF01645','PF01913','PF01993','PF02007','PF02142','PF02239','PF02240','PF02241','PF02249','PF02289','PF02315','PF02335','PF02461','PF02552','PF02560','PF02663','PF02665','PF02741','PF02745','PF02771','PF02783','PF02786','PF02787','PF02812','PF02975','PF03060','PF03063','PF03069','PF03139','PF03201','PF03404','PF03460','PF03598','PF03599','PF03951','PF04060','PF04206','PF04207','PF04208','PF04210','PF04211','PF04267','PF04268','PF04324','PF04422','PF04432','PF04879','PF04898','PF05088','PF05369','PF05440','PF06253','PF06433','PF07732','PF07969','PF07992','PF08669','PF09472','PF09505','PF09723','PF10621','PF11576','PF11844','PF12176','PF12800','PF12838','PF13187','PF13247','PF13360','PF13442','PF13447','PF13484','PF13510','PF13522','PF13806','PF14691','PF14710','PF14711')

Iron = c('PF00032','PF00033','PF00034','PF00115','PF00116','PF00355','PF00378','PF00465','PF00501','PF00510','PF00550','PF00668','PF01106','PF01154','PF01521','PF01592','PF01794','PF02167','PF02421','PF02535','PF03169','PF03264','PF07664','PF07670','PF08022','PF08030','PF08540','PF09699','PF10399','PF11854','PF13434','PF13523','PF14522')


toMatch= PfamConvert[SulfurCycle,'acc']

toMatch= PfamConvert[Oxygen,'acc']
toMatch= PfamConvert[Nitrogen,'acc']
toMatch= PfamConvert[Carbon,'acc']
toMatch= PfamConvert[Iron,'acc']
toMatch= PfamConvert[Biosynth,'acc']

toMatch= PfamConvert[Pfam2GO[grepl('cobalamin',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('flag',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('nitr',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('electron',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('aerobic',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('transport',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('fatty',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('NAD',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('biotin',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('proton',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('carbohy',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('cycle',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('tricarboxylic acid cycle',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('gluconeogenisis',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('lipid',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('virus',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('viral',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('regulation',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('trypto',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('transposase',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('oxidoreductase',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('acetyl-CoA',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('acetate',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('carbon',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('signa',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('detox',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
toMatch= PfamConvert[Pfam2GO[grepl('sugar',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']

 toMatch= PfamConvert[Pfam2GO[grepl('heme',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']

toMatch= PfamConvert[Pfam2GO[grepl('superoxi',Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']

toMatch = PfamConvert[grepl('nitr', PfamConvert$Function,ignore.case=TRUE),'acc']
toMatch = PfamConvert[grepl('NADH', PfamConvert$Function,ignore.case=TRUE),'acc']
toMatch = PfamConvert[grepl('Lactat', PfamConvert$Function,ignore.case=TRUE),'acc']
toMatch = PfamConvert[grepl('Fatty', PfamConvert$Function,ignore.case=TRUE),'acc']
toMatch = PfamConvert[grepl('Cobala', PfamConvert$Function,ignore.case=TRUE),'acc']
toMatch = PfamConvert[grepl('biotin', PfamConvert$Function,ignore.case=TRUE),'acc']
toMatch = PfamConvert[grepl('Ribosomal', PfamConvert$Function,ignore.case=TRUE),'acc']
toMatch = PfamConvert[grepl('glutathione', PfamConvert$Function,ignore.case=TRUE),'acc']
toMatch = PfamConvert[grepl('Pyruvate', PfamConvert$Function,ignore.case=TRUE),'acc']
toMatch = PfamConvert[grepl('porin', PfamConvert$Function,ignore.case=TRUE),'acc']
toMatch = PfamConvert[grepl('Flagel', PfamConvert$Function,ignore.case=TRUE),'acc']

toMatch= PfamConvert[Biosynth,'acc']


mat = t(PfamData[grepl(paste(toMatch,collapse="|"),rownames(PfamData)),])
colnames(mat)=paste(PfamConvert[PfamConvert$acc %in% colnames(mat),'ID'], colnames(mat),PfamConvert[PfamConvert$acc %in% colnames(mat),'Function'],sep=' - ')
clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mat))
combined_ordered_matrix <- mat[new_order,]

	heatmap.plus(as.matrix(t(combined_ordered_matrix)), col=colorRampPalette(c('white','black','orange','red'))(n = 150),cexRow=0.6,cexCol=0.6,scale='none',margins=c(10,20),Colv= pruned_tree_d)

toMatch = PfamConvert[grepl('transp', PfamConvert$function),'acc']


for(i in unique(Pfam2GO$GOfunc)){
	tryCatch({
	print(i)
	toMatch= PfamConvert[Pfam2GO[grepl(i,Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
	mat = t(PfamData[grepl(paste(toMatch,collapse="|"),rownames(PfamData)),])
	colnames(mat)=paste(PfamConvert[PfamConvert$acc %in% colnames(mat),'ID'], colnames(mat),PfamConvert[PfamConvert$acc %in% colnames(mat),'Function'],sep=' - ')
	clade_order <- order.dendrogram(pruned_tree_d)
	clade_name <- labels(pruned_tree_d)
	clade_position <- data.frame(clade_name,clade_order)
	clade_position <- clade_position[order(clade_position$clade_order),]
	new_order <- match(clade_position$clade_name, row.names(mat))
	combined_ordered_matrix <- mat[new_order,]
    
    pdf(paste("Pfam_per_GO_","_", i , ".pdf", sep = ""))
	heatmap.2(as.matrix(t(combined_ordered_matrix)), col=colorRampPalette(c('white','black','orange','red'))(n = 150),cexRow=0.4,cexCol=0.4,scale='none',margins=c(10,20),Colv= pruned_tree_d,trace='none',ColSideColors= as.character(cladecolorsHeat[rownames(combined_ordered_matrix)]),main=i)
	dev.off()
	for(o in c(1:10)){
		dev.off()
	}
	
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

###### need to debig this, as it is not doing what I want, the problem resides in 1 gene hits, than the thing is not converted to a matrix and rowsums does not work

SummaryGO = c()
GOtracker = c()
for(i in unique(Pfam2GO$GOfunc)){
	tryCatch({
	print(i)
	toMatch= PfamConvert[Pfam2GO[grepl(i,Pfam2GO$GOfunc,ignore.case=TRUE),'ID'],'acc']
	mat = t(PfamData[grepl(paste(toMatch,collapse="|"),rownames(PfamData)),])
	colnames(mat)=paste(PfamConvert[PfamConvert$acc %in% colnames(mat),'ID'], colnames(mat),PfamConvert[PfamConvert$acc %in% colnames(mat),'Function'],sep=' - ')
	clade_order <- order.dendrogram(pruned_tree_d)
	clade_name <- labels(pruned_tree_d)
	clade_position <- data.frame(clade_name,clade_order)
	clade_position <- clade_position[order(clade_position$clade_order),]
	new_order <- match(clade_position$clade_name, row.names(mat))
	combined_ordered_matrix <- mat[new_order,]
    
   	SummaryGO = cbind(SummaryGO, rowSums(as.matrix(mat)))
	GOtracker = c(GOtracker,i)
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

colnames(SummaryGO)= GOtracker
SummaryGO  = SummaryGO[,colSums(SummaryGO)>0]
SummaryGO  = SummaryGO[,var(SummaryGO)>0]






##########################

clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, colnames(geneCount))
combined_ordered_matrix <- geneCount[, new_order]

heatmap.plus(as.matrix(combined_ordered_matrix), RowSideColors= as.matrix(DF_significance[rownames(combined_ordered_matrix),]),col=colorRampPalette(c("white","black","black","black","black","black"))(n = 150),Colv= pruned_tree_d)

heatmap.2(as.matrix(combined_ordered_matrix),,trace="none", col=colorRampPalette(c("white","black",'orange'))(n = 150),Colv= pruned_tree_d,ColSideColors= as.character(cladecolorsHeat[colnames(combined_ordered_matrix)]),srtCol=45,cexRow=0.6,cexCol=0.6,main='')


all.sig  = unique(c(as.character(RFselect), as.character(ScoaryEnrichedDepleted$gene), as.character(Kruskal.sig.0.05)))

heatmap.plus(as.matrix(combined_ordered_matrix[all.sig,]), RowSideColors= as.matrix(DF_significance[rownames(combined_ordered_matrix[all.sig,]),]),col=colorRampPalette(c("white","black","black","black","black","black"))(n = 150),Colv= pruned_tree_d)

heatmap.2(as.matrix(combined_ordered_matrix[all.sig,]),,trace="none", col=colorRampPalette(c("white","black","black","orange",'orange','orange','red'))(n = 150),Colv= pruned_tree_d,ColSideColors= as.character(cladecolorsHeat[colnames(combined_ordered_matrix[all.sig,])]),srtCol=45,cexRow=0.6,cexCol=0.6,main='')


clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, colnames(geneCount))
combined_ordered_matrix <- geneCount[, new_order]

heatmap.plus(as.matrix(log(t(combined_ordered_matrix+1))), ColSideColors= as.matrix(DF_significance[rownames(combined_ordered_matrix[all.sig,]),]),col=colorRampPalette(c("white","black","black","black","black","black"))(n = 150),Rowv= pruned_tree_d)


heatmap.2(as.matrix(log (t(combined_ordered_matrix[all.sig,]+1))),trace="none", col=colorRampPalette(c("#30357B",'white','#DD8920','#DD8920','#DD8920'))(n = 150),Rowv= pruned_tree_d,RowSideColors= as.character(cladecolorsHeat[colnames(combined_ordered_matrix[all.sig,])]),srtCol=45,cexRow=0.6,cexCol=0.6,main='')

heatmap.plus(as.matrix(log (t(combined_ordered_matrix[all.sig,]+1))),trace="none", col=colorRampPalette(c("#30357B",'white','#DD8920','#DD8920','#DD8920'))(n = 150),Rowv= pruned_tree_d,ColSideColors= as.matrix(DF_significance[rownames(combined_ordered_matrix[all.sig,]),]),srtCol=45,cexRow=0.6,cexCol=0.6,main='')

number_of_shown_genes = as.matrix(colSums(geneCount[all.sig,]))

values2colors(number_of_shown_genes[,1], col.start="white", col.end="#442081", col.pal=NULL, na.col="grey50",breaks=10)

genome_annot_colors = data.frame(as.matrix(cladecolorsHeat))
genome_annot_colors$nr_of_genes = values2colors(number_of_shown_genes[,1], col.start="white", col.end="#442081", col.pal=NULL, na.col="grey50",breaks=10)

genome_annot_colors$GC = values2colors(meta_data$GC, col.start="white", col.end="black", col.pal=NULL, na.col="grey50")

heatmap.plus(as.matrix(log (t(combined_ordered_matrix[all.sig,]+1))),trace="none", col=colorRampPalette(c("#30357B",'white','#DD8920','#DD8920','#DD8920'))(n = 150),Rowv= pruned_tree_d,ColSideColors= as.matrix(DF_significance[rownames(combined_ordered_matrix[all.sig,]),]),srtCol=45,cexRow=0.6,cexCol=0.6,main='',RowSideColors = as.matrix(genome_annot_colors[colnames(combined_ordered_matrix[all.sig,]),]))



heatmap.2(as.matrix(log(combined_ordered_matrix+1)),trace="none"
,col=colorRampPalette(c("#30357B",'white','#DD8920','#DD8920','#DD8920'))(n = 150),srtCol=45,cexRow=0.6,cexCol=0.6,main='',margins=c(10,20),RowSideColors= cladecolors[new_order],Rowv= pruned_tree_d)




mat = matrix(rnorm(10000), nr = 1000)
rownames(mat) = sprintf("%.2f", rowMeans(mat))
subset = sample(1000, 20)


labels = rownames(mat)[subset]

Heatmap(mat, show_row_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE) + 
rowAnnotation(link = row_anno_link(at = subset, labels = labels),
  width = unit(1, "cm") + max_text_width(labels))

mat = (combined_ordered_matrix[all.sig,]+1)
subset = c("PC_00002182","PC_00003231","PC_00010167","PC_00003799","PC_00004228","PC_00004325","PC_00003832","PC_00003595","PC_00002576","PC_00006699","PC_00005046","PC_00004294","PC_00004570","PC_00004679","PC_00003780","PC_00003026","PC_00002715","PC_00002890","PC_00003588","PC_00002753")
labels = rownames(mat)[subset]

Heatmap(as.matrix(mat), show_row_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE) + 
rowAnnotation(link = row_anno_link(at = subset, labels = subset),
  width = unit(1, "cm") + max_text_width(labels))







##################
#cluster on the level of genes
#need to continue, color them accordingly and look for the cluster which is correlated with certain traits
# closer look whithin clusters to observe any peculiarities, like the algicola clade cluster, etc
# can be used to color the bars, as in figure, or to zoom into specific clusters and blaat


c1 = kmeans(combined_ordered_matrix[all.sig,], 11, nstart = 25)

#
clustercolors = values2colors((c1$cluster), col.start="green", col.end="purple", col.pal=NULL, na.col="grey50")
clustercolors_matrix =  as.matrix(c1$cluster)
clustercolors_matrix[,1] = clustercolors

gene_annot_colors = cbind(DF_significance[rownames(clustercolors_matrix),], clustercolors_matrix)



heatmap.plus(as.matrix(log (t(combined_ordered_matrix[all.sig,]+1))),trace="none", col=colorRampPalette(c("#30357B",'white','#DD8920','#DD8920','#DD8920'))(n = 150),Rowv= pruned_tree_d,ColSideColors= as.matrix(gene_annot_colors[rownames(combined_ordered_matrix[all.sig,]),]),srtCol=45,cexRow=0.6,cexCol=0.6,main='',RowSideColors = as.matrix(genome_annot_colors[colnames(combined_ordered_matrix[all.sig,]),]))

pheatmap(combined_ordered_matrix[all.sig,])


for(i in 1:11){
	heatmap.plus(as.matrix(log (t(combined_ordered_matrix[names(c1$cluster[c1$cluster==i]),]+1))),trace="none", col=colorRampPalette(c("#30357B",'white','#DD8920','#DD8920','#DD8920'))(n = 150),Rowv= pruned_tree_d,ColSideColors= as.matrix(DF_significance[rownames(combined_ordered_matrix[names(c1$cluster[c1$cluster==i]),]),]),srtCol=45,cexRow=0.6,cexCol=0.6,main='',RowSideColors = as.matrix(genome_annot_colors[colnames(combined_ordered_matrix[names(c1$cluster[c1$cluster==i]),]),]))

}

for(i in 1:11){
	heatmap.plus(as.matrix((t(combined_ordered_matrix[names(c1$cluster[c1$cluster==i]),]+1))),trace="none", col=colorRampPalette(c("#30357B",'white','#DD8920','#DD8920','#DD8920'))(n = 150),Rowv= pruned_tree_d,ColSideColors= as.matrix(DF_significance[rownames(combined_ordered_matrix[names(c1$cluster[c1$cluster==i]),]),]),srtCol=45,cexRow=0.6,cexCol=0.6,main='',RowSideColors = as.matrix(genome_annot_colors[colnames(combined_ordered_matrix[names(c1$cluster[c1$cluster==i]),]),]))

}


###########

set.seed(12345)
mat = matrix(rnorm(100), nr = 25)

hr = hclust(dist(combined_ordered_matrix[all.sig,]), method = "average")
clusters = dendextend::cutree(hr, k = 5)

Heatmap(combined_ordered_matrix[all.sig,], split = clusters,col=colorRampPalette(c("#30357B",'white','#DD8920','#DD8920','#DD8920'))(n = 150))

###########

mat = combined_ordered_matrix[all.sig,]
ha1 = HeatmapAnnotation(dist1 = anno_barplot(colSums(mat), bar_width = 1, gp = gpar(col = NA, fill = "#FFE200"), 
    border = FALSE, axis = TRUE))
ha2 = rowAnnotation(dist2 = anno_barplot(rowSums(mat), bar_width = 1, gp = gpar(col = NA, fill = "#FFE200"), 
    border = FALSE, which = "row", axis = TRUE))#, width = unit(1, "cm"))
ha_column = HeatmapAnnotation(cn = function(index) {
    year = as.numeric(colnames(mat))
    which_decade = which(year %% 10 == 0)
    grid.text(year[which_decade], which_decade/length(year), 1, just = c("center", "top"))
})
Heatmap(mat, name = "cases", col = colorRamp2(c(0, 800, 1000, 127000), c("white", "cornflowerblue", "yellow", "red")),
    cluster_columns = FALSE, show_row_dend = FALSE, rect_gp = gpar(col= "white"), show_column_names = FALSE,
    row_names_side = "left", row_names_gp = gpar(fontsize = 10),
    column_title = '',
    top_annotation = ha1, top_annotation_height = unit(1, "cm"),
    bottom_annotation = ha_column, bottom_annotation_height = grobHeight(textGrob("1900"))) + ha2

decorate_heatmap_body("cases", {2
    i = which(colnames(mat) == "1961")
    x = i/ncol(mat)
    grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 2))
    grid.text("", x, unit(1, "npc") + unit(5, "mm"))
})



##### pearson

heatmap.2(as.matrix(data.matrix(cor((combined_ordered_matrix[all.sig,]),method = 'pearson'))),trace="none",scale="none",col=colorRampPalette(c("blue4","blue4","blue4","blue4","blue4","blue4","blue4","blue4","blue4","blue4",'cyan4','cyan4','cyan4','cyan4' ,"cadetblue1","cadetblue1","cadetblue1",'darkgoldenrod1','yellow'))(n = 150),dendrogram="row", RowSideColors= as.character(cladecolorsHeat[colnames(combined_ordered_matrix[all.sig,])]))

#######################

#compare the selected dataset to the overall dataset to see whether it is representative

 bd = vegdist(t(panmatrix), method="bray")
bd_red = vegdist(t(panmatrix[all.sig,]), method="bray")

rcorr(bd, bd_red, type="pearson")$r[1,2]

jbd = vegdist(t(panmatrix), method="jaccard")
jbd_red = vegdist(t(panmatrix[all.sig,]), method="jaccard")

rcorr(jbd, jbd_red, type="pearson")$r[1,2]


fullPANTREE <- panTree(panmatrix,nboot=100)
redPANTREE <- panTree(panmatrix[all.sig,],nboot=100)




#######################



DF_signicicance = c()
for(i in rownames(t(combined_ordered_matrix))){
	DF_signicicance = rbind(DF_signicicance ,c(i , i %in% RFselect, i %in% ScoaryEnrichedDepleted$gene, i %in% Kruskal.sig.0.05))	
}
colnames(DF_signicicance) = c('gene','RF','Scoary','KruskalFDR')

for(i in ANVIOfiltered$protein_cluster_id){
	if(i %in% DF_signicicance$gene){}
	else{
		rbind(DF_signicicance ,c(i , i %in% RFselect, i %in% ScoaryEnrichedDepleted$gene, i %in% Kruskal.sig.0.05))
		}
}

DF_signicicance[DF_signicicance=='TRUE'] = 'black'
DF_signicicance[DF_signicicance=='FALSE'] = 'gray'
rownames(DF_signicicance) = DF_significance$gene
DF_signicicance = DF_signicicance[,2:length(colnames(DF_signicicance))]
heatmap.plus(as.matrix(t(combined_ordered_matrix)), RowSideColors= DF_signicicance)



rownames(t(combined_ordered_matrix))


###########


intersect = intersect(Kruskal.sig, RFselect)

bd = vegdist(t(panmatrix), method="bray") 
Jd = distJaccard(t(panmatrix))

venn(list('KruskalFDR' = kruskal.wallis.table[kruskal.wallis.table$q.value < 0.05 ,c('id')], 'Random Forest' = RFselect))
ScoaryEnrichedDepleted

venn(list('KruskalFDR' = kruskal.wallis.table[kruskal.wallis.table$q.value < 0.05 ,c('id')], 'Random Forest' = RFselect,'Scoary'= ScoaryEnrichedDepleted$gene))
intersect = intersect(Kruskal.sig, RFselect,ScoaryEnrichedDepleted$gene)

venn(list('KruskalFDR' = Kruskal.sig.0.05, 'Random Forest' = RFselect,'Scoary'= ScoaryEnrichedDepleted$gene))
intersect = intersect(Kruskal.sig.0.05, as.character(RFselect), as.character(ScoaryEnrichedDepleted$gene))


 intersect = intersect( as.character(RFselect), as.character(ScoaryEnrichedDepleted$gene))
intersect = intersect(Kruskal.sig.0.05, intersect)

Kruskal.sig = kruskal.wallis.table[kruskal.wallis.table$q.value < 0.05 ,c('id')]

intersect_bd = vegdist(t(panmatrix[intersect,]), method="bray") 
intersect_Jd = distJaccard(t(panmatrix[intersect,]))

intersect2_bd = vegdist(t(panmatrix[intersect(kruskal.wallis.table[kruskal.wallis.table$q.value < 0.05 ,c('id')], RFselect),]), method="bray") 
intersect2_Jd = distJaccard(t(panmatrix[intersect(kruskal.wallis.table[kruskal.wallis.table$q.value < 0.05 ,c('id')], RFselect),]))


RF_bd = vegdist(t(panmatrix[RFselect,]), method="bray") 
RF_Jd = distJaccard(t(panmatrix[RFselect,]))


rcorr(bd, intersect_bd, type="pearson")$r[1,2]
rcorr(Jd, intersect_Jd, type="pearson")$r[1,2]

rcorr(bd, intersect2_bd, type="pearson")$r[1,2]
rcorr(Jd, intersect2_Jd, type="pearson")$r[1,2]

rcorr(bd, RF_bd, type="pearson")$r[1,2]
rcorr(Jd, RF_Jd, type="pearson")$r[1,2]



panmatrix


bd = vegdist(t(panmatrix), method="bray") 
Jd = distJaccard(t(panmatrix))

mbd <- melt(as.matrix(bd))[melt(upper.tri(as.matrix(bd)))$value,]
mJd <- melt(as.matrix(Jd))[melt(upper.tri(as.matrix(Jd)))$value,]
rcorr(cbind(mbd$value, mJd$value), type="pearson")
rcorr(cbind(mbd$value, mJd$value), type="pearson")$r[1,2]

bd = vegdist(t(panmatrix), method="bray") 

reduced_bd = vegdist(t(panmatrix[1:1000,]), method="bray") 
reduced_bd_r <- melt(as.matrix(reduced_bd))[melt(upper.tri(as.matrix(reduced_bd)))$value,]
rcorr(cbind(mbd$value, reduced_bd_r$value), type="pearson")$r[1,2]

ckr = c()
for(i in 2:nrow(panmatrix)){
	print(i)
	reduced_bd = vegdist(t(panmatrix[1: i,]), method="bray") 
	reduced_bd_r <- melt(as.matrix(reduced_bd))[melt(upper.tri(as.matrix(reduced_bd)))$value,]
	ckr = rbind(ckr ,c(rcorr(cbind(mbd$value, reduced_bd_r$value), type="pearson")$r[1,2]))

}
plot(ckr)

plot(db, Jd)




#include one RF-determined variable at the time and compare the distance matrices to the original dataset
incremental_set = c(var_importance[1,]$variable)
incr = c()
for(i in 2:nrow(var_importance)){
	incremental_set = c(incremental_set,var_importance[i,]$variable)
	print(i)
	reduced_bd = vegdist(t(panmatrix[incremental_set,]), method="bray") 
	incr = rbind(incr ,c(rcorr(reduced_bd, bd, type="pearson")$r[1,2]))

}
plot(incr)

plot(db, Jd)



###################################################
# using phyloseq package, has function to plot microbiome networks
###################################################

library(phyloseq); packageVersion("phyloseq")


panmatrix2 = as.matrix(geneCount)

panmatrix2 = as.matrix(panmatrix)
(panmatrix2 <- unclass(panmatrix2))

OTU = otu_table(panmatrix2, taxa_are_rows = TRUE)
OTUpa = OTU
OTUpa[OTUpa>0] <-1

OTUpa = otu_table(OTUpa, taxa_are_rows = TRUE)


sample_names()

sampledata = sample_data(grouping_factors)
sample_names(sampledata)=
phyloseq = phyloseq(OTU, sampledata)
phyloseq2 = phyloseq(OTUpa, sampledata)

set.seed(711L)


ig <- make_network(phyloseq,dist.fun="jaccard")
plot_network(ig, phyloseq, color="group2", line_weight=1, label=NULL,hjust=1)


ig <- make_network(phyloseq2,distance="jaccard")
plot_network(ig, phyloseq2, color="group2", line_weight=1, label=NULL,hjust=1)

ig <- make_network(phyloseq2,distance="bray")
plot_network(ig, phyloseq2, color="group2", line_weight=1, label=NULL,hjust=1)


ig <- make_network(phyloseq2,distance="jaccard",type='taxa')
plot_network(ig,type='taxa', phyloseq2,  line_weight=1, label=NULL,hjust=1)




dist_methods = dist_methods[-which(dist_methods=="ANY")]
dist_methods <- unlist(distanceMethodList)

dist_methods = dist_methods[c(5,6,8,10,11,17)]
#dist_methods[(1:2)]


plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
	print(paste('calculating ordination with distance method:',i))
    # Calculate distance matrix
    iDist <- distance(phyloseq, method=i)
    # Calculate ordination
    #iMDS  <- ordinate(phyloseq, "MDS", distance=iDist)
    iMDS  <- ordinate(phyloseq, "MDS", distance=iDist)
    
    ## Make plot
    # Don't carry over previous plot (if error, p will be blank)
    p <- NULL
    # Create plot, store as temp variable, p
    p <- plot_ordination(phyloseq, iMDS, color="group2")
    # Add title to each plot
    p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
    # Save the graphic to file.
    plist[[i]] = p
}


df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=group2))
p = p + geom_point(size=2) + scale_colour_manual(values = c("#39811D","#86DA81","#A6D8D4","#F2A968","#F2EC70","#E38FDD","#898989","#76AECF","#B34D22"))
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics")
p



plot_richness(phyloseq, x="group2", color ="group2") + scale_colour_manual(values = c("#39811D","#86DA81","#A6D8D4","#F2A968","#F2EC70","#E38FDD","#898989","#76AECF","#B34D22"))+geom_boxplot()



##########################################################################################################################
#
#
#	BIOLOGICAL MARKERS - PROOF OF CONCEPT
#
#
##########################################################################################################################
# the general idea is here to test the biomarkers on a wider set of genomes 

##########################


mat = t(table(droplevels(subset(ANVIO, protein_cluster_id  %in% all.sig)[,c('protein_cluster_id','genome_name')])))

ANVIO_cat2$quality = ifelse(ANVIO_cat2$Bin_Id %in% colnames(geneCount),"included", "mapped")
ANVIO_cat2$qualityColor = ifelse(ANVIO_cat2$Bin_Id %in% colnames(geneCount),"black", "red")


heatmap.2(as.matrix(mat),trace='none',col=colorRampPalette(c('white',"#000033", "#FF3300", "#FF3300", "#FF3300"))(n = 50),margins=c(12,12),RowSideColors=ANVIO_cat2$qualityColor)

mat <- unclass(mat)

fullcount.sig = as.matrix(t(mat))
#panmatrix2 = as.matrix(panmatrix)
#(panmatrix2 <- unclass(panmatrix2))

OTU = otu_table(fullcount.sig ,taxa_are_rows = TRUE)
#OTUpa = OTU
#OTUpa[OTUpa>0] <-1

#OTUpa = otu_table(OTUpa, taxa_are_rows = TRUE)
#grouping_factors2 = ANVIO_cat2
#grouping_factors2 = grouping_factors[grouping_factors2$Bin_Id %in% colnames(fullcount.sig),]
#rownames(grouping_factors2)= grouping_factors2$Bin_Id
#grouping_factors2 = grouping_factors2[,2:length(colnames(grouping_factors2))]

sampledata = sample_data(ANVIO_cat2)

sample_names()
sample_names(sampledata) = sample_names(OTU)

phyloseq = phyloseq(OTU, sampledata)

    iMDS  <- ordinate(phyloseq, "MDS", distance='bray')
    p <- plot_ordination(phyloseq, iMDS, color="group2",shape='quality')
    # Add title to each plot
    p <- p + ggtitle(paste("MDS - bray distance method", sep=""))+ geom_point(size=3) + scale_colour_manual(values = colors)

    iMDS  <- ordinate(phyloseq, "MDS", distance='manhattan')
    pj <- plot_ordination(phyloseq, iMDS, color="group2",shape='quality')
    # Add title to each plot
    pj <- pj + ggtitle(paste("MDS - manhattan distance method", sep=""))+ geom_point(size=3) + scale_colour_manual(values = colors)

multiplot(p,pj,cols=2)


##### DISTANCE BASED
plot(dist(t(fullcount.sig)))
hclustdistfull.sig = hclust(vegdist(t(fullcount.sig), method="bray"), method = "average")
dist.full.sug = vegdist(t(fullcount.sig), method="bray")
dist.full.sug.jaccard = vegdist(t(fullcount.sig), method="jaccard")

dndfull.sig <- as.dendrogram(hclustdistfull.sig)
meltfull.sig <- melt(as.matrix(dist.full.sug))[melt(upper.tri(as.matrix(dist.full.sug)))$value,]
meltfull.sig.jaccrard <- melt(as.matrix(dist.full.sug.jaccard))[melt(upper.tri(as.matrix(dist.full.sug.jaccard)))$value,]

tanglegram(rep_tree3_d, dndfull.sig)

PhyloDistMatrix3<-cophenetic(tree3)
PhyloDistMatrix3  = PhyloDistMatrix3[ order(row.names(PhyloDistMatrix3)), ]
PhyloDistMatrix3  = PhyloDistMatrix3[ , order(colnames(PhyloDistMatrix3))]
melt4 <- melt(as.matrix(PhyloDistMatrix3))[melt(upper.tri(as.matrix(PhyloDistMatrix3)))$value,]

comparePhylo = data.frame('RAxMLCampbellRooted' = melt4$value, 'functional.bray' = meltfull.sig$value,'functional.jaccard' = meltfull.sig.jaccrard $value)
bp = ggscatter(comparePhylo, x = "RAxMLCampbellRooted", y = "functional.bray", add = "loess", conf.int = TRUE,main='Bray')
jp = ggscatter(comparePhylo, x = "RAxMLCampbellRooted", y = "functional.jaccard", add = "loess", conf.int = TRUE,main='Jaccard')
multiplot(bp,jp,cols=2)


##########################################################################################################################

Mapptree <- panTree(t(fullcount.sig),nboot=1000)
fullCampbel = read.tree("~/DATA/MarinobacterGenomics/2018_ProkComp/RAxML_bipartitions.Campbell_AA")


ladderize(RAxMLANVIORooted, right = FALSE)

tr1 = ggtree(midpoint.root(as.phylo(Mapptree $Htree)),layout='circular') %<+% ANVIO_cat2 + geom_point(aes(color = group2,shape=quality,size=3)) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2)+ggtitle('protein clusters')+geom_treescale()
tr2 = ggtree(ladderize(midpoint.root(fullCampbel), right = FALSE),branch.length='none',layout='circular')%<+% ANVIO_cat2 + geom_point(aes(color = group2,shape=quality,size=3)) + scale_colour_manual(values=colors)+ggtitle('phylogeny')+geom_treescale()


multiplot(tr1,tr2,cols=2)

tr1 = open_tree(tr1, angle=180)
tr2 = open_tree(tr2, angle=180)
multiplot(tr1,tr2,cols=2)

tr1 = ggtree(midpoint.root(as.phylo(Mapptree $Htree)),layout='circular') %<+% ANVIO_cat2 + geom_point(aes(color = group2,shape=quality,size=2)) + scale_colour_manual(values=colors)+ geom_text2(aes(label=label, subset=!isTip), hjust=1.2,vjust=-.5,size=2)+ggtitle('protein clusters')+geom_treescale()
tr1 = open_tree(tr1, angle=180)


ANVIOmarkers = subset(ANVIO, protein_cluster_id  %in% all.sig)

Markers_algicola = ANVIOmarkers[ANVIOmarkers$genome_name == 'Marinobacter_algicola_DG893',]






##########################################################################################################################
#
#	Negative binomial GLM fitting and Wald statistics for abundance data using DESeq {DESeq2} package
#
#		* Differential expression analysis based on negative binomial (a.k.a gamma-distribution)
#		* Tailored to identify group specific genes from a pangenome-matrix (any matrix)(nxp)
#
##########################################################################################################################

library(DESeq2)

#convert the table to DESeqDataSet object
geneCount = round(as(panmatrix, "matrix"), digits = 0)

# We will add 1 to the countData otherwise DESeq will fail with the error:
# 	estimating size factors
# 	Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  :
# 	every gene contains at least one zero, cannot compute log geometric means

geneCount <- geneCount +1 
rownames(grouping_factors)=NULL

dds <- DESeqDataSetFromMatrix(geneCount, grouping_factors, as.formula(~ lifestyle))
data_deseq_test = DESeq(dds)

## Extract the results
res = results(data_deseq_test, cooksCutoff = FALSE)
res_tax = cbind(as.data.frame(res), as.matrix(geneCount[rownames(res), ]), OTU = rownames(res))

sig = 0.05
fold = 0
plot.point.size = 2
label=F
tax.display = NULL
tax.aggregate = "OTU"

res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange))
res_tax_sig <- res_tax_sig[order(res_tax_sig$padj),]

res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
res_tax$Significant[is.na(res_tax$Significant)] <- "No"
p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) +
 geom_point(size = plot.point.size) +
 scale_x_log10() +
 scale_color_manual(values=c("black", "red")) +
 labs(x = "Mean abundance", y = "Log2 fold change")+theme_bw()
if(label == T){
 if (!is.null(tax.display)){
   rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
 } else {
   rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
 }
 p1 <- p1 + geom_text(data = subset(rlab, Significant == "Yes"), aes(label = Display), size = 4, vjust = 1)
}
pdf("NB_MA.pdf")
print(p1)
dev.off()

res_tax_sig_abund = cbind(as.data.frame(countData[rownames(res_tax_sig), ]), OTU = rownames(res_tax_sig), padj = res_tax[rownames(res_tax_sig),"padj"]) 

#Apply normalisation (either use relative or log-relative transformation)
#data<-abund_table/rowSums(abund_table)
data<-log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
data<-as.data.frame(data)

#Now we plot taxa significantly different between the categories
df<-NULL
for(i in res_tax[rownames(res_tax_sig),"OTU"]){
 tmp<-data.frame(data[,i],grouping_info$Country,rep(paste(paste(i,gsub(".*;","",gsub(";+$","",paste(sapply(OTU_taxonomy[i,],as.character),collapse=";"))))," padj = ",sprintf("%.5g",res_tax[i,"padj"]),sep=""),dim(data)[1]))
 if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Type","Taxa")

p<-ggplot(df,aes(Type,Value,colour=Type))+ylab("Log-relative normalised")
p<-p+geom_boxplot()+geom_jitter()+theme_bw()+
 facet_wrap( ~ Taxa , scales="free_x",nrow=1)
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90))
pdf("NB_significant.pdf",width=160,height=10)
print(p)



fvg = vegdist(t(geneCount), method="bray")
rvg = vegdist(t(geneCount[all.sig,]), method="bray")
plot(fvg,rvg)


dd <- as.dist((cor(((geneCount)), method='spearman')))
dd= as.matrix(dd)
dd[dd==0]=1

heatmap.2(dd,col = (colorRampPalette(brewer.pal(9, "BrBG"))(100)),trace='none',RowSideColors= cladecolorsHeat)



ff <- as.dist((cor(((geneCount[all.sig,])), method='spearman')))
ff= as.matrix(ff)
ff[ff==0]=1

heatmap.2(ff,col = (colorRampPalette(brewer.pal(9, "BrBG"))(100)),trace='none',RowSideColors= cladecolorsHeat)


gg <- as.dist((cor(((geneCount[intersect,])), method='spearman')))
gg= as.matrix(gg)
gg[gg==0]=1

heatmap.2(gg,col = (colorRampPalette(brewer.pal(9, "BrBG"))(100)),trace='none',RowSideColors= cladecolorsHeat)
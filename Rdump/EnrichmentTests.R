# 2018 Frederik De Boever


################################
#		Dependencies
################################
library(ape)
library(phytools)
library(phylolm)
library(gplots)

print ("algicola versus hydrocarbonoclasticus")
tree<-read.tree("~/DATA/MarinobacterGenomics/2018_ProkComp/SCO_682_RAxML_bipartitions_autosubst_b100")
tree<- midpoint.root(tree)

dummy = data.frame(genome = levels(ANVIO$genome))

Map<-ANVIO_cat
Tab<-AnvioData
Tab = t(Tab)
write.csv(Tab , file= paste('~/Tab','_added.results.csv',sep=''),row.names = FALSE)
write.csv(Map , file= paste('~/Map','_added.results.csv',sep=''),row.names = FALSE)


################################
#		RUN ANALYSIS
################################

#hypergeometric test
finRes2 = HyperGeometric_test(tree,Tab,Map,'algicola')
GroupA <-finRes2[which(finRes2$p.value_enriched_binary<=0.05),"protein_cluster_id"]

#phyloGLM
df_traits = Select_trait_phyloGLM(Map,'algicola')
Prepare_data_phyloGLM(Tab, df_traits,'algicola')
Create_Dir("TESTRUN","split_matrix_pa_npa_")
perform_PhyloGLM(read.table('~/split_matrix_pa_npa_TESTRUN/splitted_file_1.tsv',header=T,sep="\t",row.names=1,check.names=F), subtree)



################################################################################################################################
#	
#					Hypergeometric test Functions
#	
################################################################################################################################

HyperGeometric_test <- function(tree,Tab,Map,group){
 	print(group)
 	Map<-droplevels(Map)
	Tab<-Tab[match(Map$name,rownames(Tab)),]
	Tab<-Tab[,which(colSums(Tab)!=0)]
 	subtree<-drop.tip(phy=tree,tip=which(!(tree$tip.label%in%rownames(Tab))))

 	finRes<-allgenes_hyp(subtree,Tab,Map, group)
    colnames(finRes)[1] = "protein_cluster_id"
    finRes2 = right_join(finRes, DEDUPLICATED[,c('protein_cluster_id','COG_CATEGORY_ACC','COG_CATEGORY','COG_FUNCTION_ACC','COG_FUNCTION')],by=('protein_cluster_id'))
 	finRes2 =  finRes2[order(-finRes2$score_enriched_binary),]
	return(finRes2)
}

allgenes_hyp <- function(tree,genes,Map,SelectedGroup){  
  N.pa <- sum(genes[ Map$group == SelectedGroup, ])
  N.npa <- sum(genes[ Map$group != SelectedGroup, ])
  #N.npa <- sum(genes[ Map$group == "algicola", ])

  Res <- NULL
  for(gene in 1:ncol(genes)){
    #gene <- 1
    gene.id <- colnames(genes)[gene]
    Map$gene <- genes[,gene]
    
 
    # Binary version
    pval <- phyper(q = nrow(subset(Map, group == SelectedGroup & gene > 0)) - 1,
                   m = sum(Map$gene > 0),
                   n = nrow(subset(Map,gene == 0)),
                   k = nrow(subset(Map, group == SelectedGroup)),lower.tail=FALSE)
    score <- -log10(pval)
    
    pval2 <- phyper(q = sum(subset(Map, group == SelectedGroup)$gene) - 1, 
                    m = N.pa, n = N.npa, k = sum(Map$gene),lower.tail=FALSE)
    
    #Test for depletion binary version 
    pval_dep <- phyper(q = nrow(subset(Map, group == SelectedGroup & gene > 0)),
                   m = sum(Map$gene > 0),
                   n = nrow(subset(Map,gene == 0)),
                   k = nrow(subset(Map, group == SelectedGroup)),
                   lower.tail=TRUE)
    score_dep<--log10(pval_dep)

    #Test for depletion Raw counts version 
    pval2_dep <- phyper(q = sum(subset(Map, group == SelectedGroup)$gene),
                    m = N.pa, n = N.npa, k = sum(Map$gene),lower.tail=TRUE)
    score2_dep<--log10(pval2_dep)

    res <- data.frame(gene.id = gene.id, score_enriched_binary = score, p.value_enriched_binary = pval,
                      score_depletion_binary=score_dep,p.value_depletion_binary=pval_dep,
                      score_enriched_rawcounts = -log10(pval2), p.value_enriched_rawcounts = pval2,
                      score_depletion_rawcounts=score2_dep,p.value_depletion_rawcounts=pval2_dep)

    Map$gene <- NULL
    Res <- rbind(Res,res)
  }
  
  Res<-data.frame(gene.id=Res$gene.id,score_enriched_binary=Res$score_enriched_binary,
  z.score_enriched_binary= (Res$score_enriched_binary - mean(Res$score_enriched_binary)) / sd(Res$score_enriched_binary),
  p.value_enriched_binary=Res$p.value_enriched_binary,
  score_depletion_binary=Res$score_depletion_binary,
  z.score_depletion_binary=(Res$score_depletion_binary - mean(Res$score_depletion_binary)) / sd(Res$score_depletion_binary),
  p.value_depletion_binary=Res$p.value_depletion_binary,
  score_enriched_rawcounts=Res$score_enriched_rawcounts,
  z.score_enriched_rawcounts=(Res$score_enriched_rawcounts - mean(Res$score_enriched_rawcounts)) / sd(Res$score_enriched_rawcounts),
  p.value_enriched_rawcounts=Res$p.value_enriched_rawcounts,
  score_depletion_rawcounts=Res$score_depletion_rawcounts,
  z.score_depletion_rawcounts=(Res$score_depletion_rawcounts - mean(Res$score_depletion_rawcounts)) / sd(Res$score_depletion_rawcounts),
  p.value_depletion_rawcounts=Res$p.value_depletion_rawcounts)

  return(Res)
}

################################################################################################################################
#	
#	PHYLOGLM related functions
#	
################################################################################################################################

Create_Dir <- function(prefix,sufix){
	outdir<-paste(prefix, sufix,"/",sep="")
	dir.create(outdir)
	return(outdir)
}

Select_trait_phyloGLM <- function(GroupingTable,group){
	tratos = data.frame( TraitY = ifelse(GroupingTable$group==group, 1, 0))
	rownames(tratos)<-GroupingTable$name
	df_traits = tratos
	return(df_traits)
}

Prepare_data_phyloGLM <- function(AbundanceMatrix, df_traits,group,outdir){
	numsplits<-as.numeric(1000000000000)
	span=numsplits-1
	num_paral=0;
	start <- seq(from=1,to=ncol(AbundanceMatrix), by = numsplits)

	for(index in start){
		num_paral=num_paral+1
		outfile=paste(outdir,"splitted_file_",num_paral,".tsv",sep="")
		top<-span + index
		if(top >ncol(AbundanceMatrix)){
			top<-ncol(AbundanceMatrix)
			subTab<-AbundanceMatrix[,index:top]
			subTab<-cbind(df_traits,as.data.frame.matrix(subTab))
			write.table(file=outfile,x=subTab,quote=F,row.names=T,col.names=T,sep="\t",append=F)
		}else{
			subTab<-AbundanceMatrix[,index:top]
			subTab<-cbind(df_traits,subTab)
			write.table(file=outfile,x=subTab,quote=F,row.names=T,col.names=T,sep="\t",append=F)
		}
}}


SplitTab<-read.table('~/split_matrix_pa_npa_TESTRUN/splitted_file_1.tsv',header=T,sep="\t",row.names=1,check.names=F)
#SplitTab = Tab
tree = midRoot_raxml_ANVIO_SCO_AA
subtree = tree

perform_PhyloGLM <- function(SplitTab, subtree){
	start=2
	y<-SplitTab[,1]
	names(y)<-rownames(SplitTab)
	Res <- NULL

	for(column in start:ncol(SplitTab)){
		print(column)
		x<-SplitTab[,column]
		orthogroup<-colnames(SplitTab)[column]
		names(x)<-rownames(SplitTab)
		dat<-as.data.frame(cbind(y,x))
	
		m1 <- tryCatch(phyloglm(formula=y~x, data = dat,phy = subtree,method = "logistic_IG10"), error = function(e) list(coefficients = NA))
    	if(is.na(coef(m1)[1])){
      	res <- data.frame(orthogroup.id = orthogroup,
                    	    Estimate = NA,
                	        SE = NA,
            	            z.value = NA,
        	                p.value = NA)
    	}else{
      	m1.sum <- summary(m1)
      	res <- data.frame(orthogroup.id = orthogroup,
                    	    Estimate = m1.sum$coefficients["x",1],
                	        SE = m1.sum$coefficients["x",2],
            	            z.value = m1.sum$coefficients["x",3],
        	                p.value = m1.sum$coefficients["x",4])
    	}
    	rm(m1,m1.sum)
    	Res <- rbind(Res,res)
	}
	return(Res)
}

outfile<-paste("~/split_matrix_pa_npa_TESTRUN/splitted_file_1",".phyloglm.tab",sep="")
write.table(x=Res,file=outfile,col.names=T,row.names=F,quote=F,append=F,sep="\t")


 
GroupA <-finRes2[which(finRes2$p.value_enriched_binary<=0.05),"protein_cluster_id"]

Map<-ANVIO_cat
Tab<-AnvioData


################################################################################################################################
#	
#	select a subtree from the data
#	
################################################################################################################################

Map<-Map[which(Map$group %in% c("algicola","hydrocarbo")),]
Map<-droplevels(Map)
Tab<-Tab[match(Map$name,rownames(Tab)),]
Tab<-Tab[,which(colSums(Tab)!=0)]
subtree<-drop.tip(phy=tree,tip=which(!(tree$tip.label%in%rownames(Tab))))







finRes<-allgenes_hyp(subtree,Tab,Map)
 colnames(finRes)[1] = "protein_cluster_id"
 
 
finRes2 = right_join(finRes, DEDUPLICATED[,c('protein_cluster_id','COG_CATEGORY_ACC','COG_CATEGORY','COG_FUNCTION_ACC','COG_FUNCTION')],by=('protein_cluster_id'))
 
finRes2 =  finRes2[order(-finRes2$score_enriched_binary),]
 
write.table(finRes2,file = "hyp_res_pa_npa.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE,append=FALSE)

COGselect = head(finRes2$COG_FUNCTION_ACC,100)
COGselect = as.character(COGselect)
COGselect = COGselect[COGselect!='']
mat = t(table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% COGselect)[,c('COG_FUNCTION','genome_name')])))














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

######PA versus NPA######
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Map<-ANVIO_cat
Tab<-AnvioData
df_trait=cbind()
for(group in unique(ANVIO_cat$group)){
	trait<-as.numeric(ifelse(Map$group==group,1,0))
	df_trait= cbind(df_trait,trait)
}
rownames(df_trait)<-Map$name
colnames(df_trait)<-unique(ANVIO_cat$group)


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

#Abundance_Table
df_Tab<-cbind(DF4,as.data.frame.matrix(Tab))

#Presence_Absence Table
Tab1 = as.data.frame.matrix(Tab)
Tab1[Tab1>0] <- 1
df_Tab_pres_abs<-cbind(DF4,Tab1)

write.table(df_trait,"traits_scoary.tsv",row.names=T,col.names=NA,sep=",",append=F,quote=F)
write.table(df_Tab,"matrix_scoary_abundance.tsv",row.names=F,col.names=T,sep=",",append=F,quote=F)
write.table(df_Tab_pres_abs,"matrix_scoary_Presence_Absence.tsv",row.names=F,col.names=T,sep=",",append=F,quote=F)

# /Users/sa01fd/Genomics/miniconda2/bin/scoary -t /Users/sa01fd/traits_scoary.tsv -g /Users/sa01fd/matrix_scoary.tsv -n /Users/sa01fd/DATA/MarinobacterGenomics/2018_ProkComp/SCO_682_RAxML_bipartitions_autosubst_b100 -e 1000 -p 1.0 --threads 8 -s 6 -o ./TEST4_SCOARY

##############################################
#	COG CATEGORY LIST
##############################################

COG_meaning = data.frame(rbind(c('A','RNA processing and modification'),c('B','Chromatin Structure and dynamics'),c('C','Energy production and conversion'),c('D','Cell cycle control and mitosis'),c('E','Amino Acid metabolis and transport'),c('F','Nucleotide metabolism and transport'),c('G','Carbohydrate metabolism and transport'),c('H','Coenzyme metabolis'),c('I','Lipid metabolism'),c('J','Tranlsation'),c('K','Transcription'),c('L','Replication and repair'),c('M','Cell wall/membrane/envelop biogenesis'),c('N','Cell motility'),c('O','Post-translational modification, protein turnover, chaperone functions'),c('P','Inorganic ion transport and metabolism'),c('Q','Secondary Structure'),c('T','Signal Transduction'),c('U','Intracellular trafficing and secretion'),c('Y','Nuclear structure'),c('Z','Cytoskeleton'),c('R','General Functional Prediction only'),c('S','Function Unknown'),c('V','Defence Mechanisms'),c('W','Extracellular structures'),c('Y','Nuclear structure'),c('X','Mobilome: prophages, transposons')))
names(COG_meaning)=c("Var1","long")



##############################################
#	SCOARY OUTPUT ANALYSER
##############################################

df_CpcogFull = cbind()
COG_distribution_all = data.frame(cbind(Var1=as.character(unique(ANVIO$COG_CATEGORY_ACC))))
for(group in unique(ANVIO_cat$group)){
	scoary <- read.table(paste('~/v4_SCOARY/',group,'_02_03_2018_1222.results.csv',sep=''),header=TRUE,sep=",") #protein clusters
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
}

names(df_CpcogFull)= c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
rownames(df_CpcogFull)=as.character(unique(ANVIO_cat$group))
df_CpcogFull = data.frame(t(df_CpcogFull))
df_CpcogFull$Var1 = c( "C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )

df_CpcogFull = right_join(COG_meaning, df_CpcogFull,by=('Var1'))
df_CpcogFull2 = df_CpcogFull[,3:ncol(df_CpcogFull)]
rownames(df_CpcogFull2) = paste(df_CpcogFull$Var1 , df_CpcogFull$long, sep=' - ')

heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=viridis(50),margins=c(12,12), srtCol =45)
heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=colorRampPalette(c("#F7FBFF", "#08306B"))(n = 50),margins=c(12,12), srtCol =45 ,cexRow = 0.8, cexCol = 0.8)

GroupB <-as.character(scoary$genes)




#######################
#	ANALYLISIS DEPLETED
#######################

df_CpcogFull = cbind()
COG_distribution_all = data.frame(cbind(Var1=as.character(unique(ANVIO$COG_CATEGORY_ACC))))
for(group in unique(ANVIO_cat$group)){
	scoary <- read.table(paste('~/v4_SCOARY/',group,'_02_03_2018_1222.results.csv',sep=''),header=TRUE,sep=",") #protein clusters
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
}

names(df_CpcogFull)= c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
rownames(df_CpcogFull)=as.character(unique(ANVIO_cat$group))
df_CpcogFull = data.frame(t(df_CpcogFull))
df_CpcogFull$Var1 = c( "C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )

df_CpcogFull = right_join(COG_meaning, df_CpcogFull,by=('Var1'))
df_CpcogFull2 = df_CpcogFull[,3:ncol(df_CpcogFull)]
rownames(df_CpcogFull2) = paste(df_CpcogFull$Var1 , df_CpcogFull$long, sep=' - ')

heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=viridis(50),margins=c(12,12), srtCol =45)
heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=colorRampPalette(c("#FFF5F0", "#67000D"))(n = 50),margins=c(12,12), srtCol =45,cexRow = 0.8, cexCol = 0.8)

#######################
#	SENSITIVITY == 100
#######################

df_CpcogFull = cbind()
COG_distribution_all = data.frame(cbind(Var1=as.character(unique(ANVIO$COG_CATEGORY_ACC))))
for(group in unique(ANVIO_cat$group)){
	scoary <- read.table(paste('~/v4_SCOARY/',group,'_02_03_2018_1222.results.csv',sep=''),header=TRUE,sep=",") #protein clusters
	scoary = right_join(DEDUPLICATED3, scoary,by=('genes'))
	scoary = scoary[order(scoary$Benjamini_H_p),]
	#scoary = scoary[scoary$Benjamini_H_p <= 0.01,]
	
	#Crucual line to change
	scoary = scoary[scoary$Sensitivity == 100,]

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
}

names(df_CpcogFull)= c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
rownames(df_CpcogFull)=as.character(unique(ANVIO_cat$group))
df_CpcogFull = data.frame(t(df_CpcogFull))
df_CpcogFull$Var1 = c( "C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )

df_CpcogFull = right_join(COG_meaning, df_CpcogFull,by=('Var1'))
df_CpcogFull2 = df_CpcogFull[,3:ncol(df_CpcogFull)]
rownames(df_CpcogFull2) = paste(df_CpcogFull$Var1 , df_CpcogFull$long, sep=' - ')

heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=viridis(50),margins=c(12,12), srtCol =45)
heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=colorRampPalette(c("#F7FBFF", "#08306B"))(n = 50),margins=c(12,12), srtCol =45 ,cexRow = 0.8, cexCol = 0.8)

#######################
#	SENSITIVITY== 0
#######################

df_CpcogFull = cbind()
COG_distribution_all = data.frame(cbind(Var1=as.character(unique(ANVIO$COG_CATEGORY_ACC))))
for(group in unique(ANVIO_cat$group)){
	scoary <- read.table(paste('~/v4_SCOARY/',group,'_02_03_2018_1222.results.csv',sep=''),header=TRUE,sep=",") #protein clusters
	scoary = right_join(DEDUPLICATED3, scoary,by=('genes'))
	scoary = scoary[order(scoary$Benjamini_H_p),]
	#scoary = scoary[scoary$Benjamini_H_p <= 0.01,]
	
	#Crucual line to change
	scoary = scoary[scoary$Sensitivity == 0,]

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
}

names(df_CpcogFull)= c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
rownames(df_CpcogFull)=as.character(unique(ANVIO_cat$group))
df_CpcogFull = data.frame(t(df_CpcogFull))
df_CpcogFull$Var1 = c( "C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )

df_CpcogFull = right_join(COG_meaning, df_CpcogFull,by=('Var1'))
df_CpcogFull2 = df_CpcogFull[,3:ncol(df_CpcogFull)]
rownames(df_CpcogFull2) = paste(df_CpcogFull$Var1 , df_CpcogFull$long, sep=' - ')

heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=viridis(50),margins=c(12,12), srtCol =45)
heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=colorRampPalette(c("#FFF5F0", "#67000D"))(n = 50),margins=c(12,12), srtCol =45,cexRow = 0.8, cexCol = 0.8)


#######################

#######################
#	NON-stringent ENRICHED
#######################

df_CpcogFull = cbind()
COG_distribution_all = data.frame(cbind(Var1=as.character(unique(ANVIO$COG_CATEGORY_ACC))))
for(group in unique(ANVIO_cat$group)){
	scoary <- read.table(paste('~/v4_SCOARY/',group,'_02_03_2018_1222.results.csv',sep=''),header=TRUE,sep=",") #protein clusters
	scoary = right_join(DEDUPLICATED3, scoary,by=('genes'))
	scoary = scoary[order(scoary$Benjamini_H_p),]
	scoary = scoary[scoary$Naive_p <= 0.01,]
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
}

names(df_CpcogFull)= c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
rownames(df_CpcogFull)=as.character(unique(ANVIO_cat$group))
df_CpcogFull = data.frame(t(df_CpcogFull))
df_CpcogFull$Var1 = c( "C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )

df_CpcogFull = right_join(COG_meaning, df_CpcogFull,by=('Var1'))
df_CpcogFull2 = df_CpcogFull[,3:ncol(df_CpcogFull)]
rownames(df_CpcogFull2) = paste(df_CpcogFull$Var1 , df_CpcogFull$long, sep=' - ')

heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=viridis(50),margins=c(12,12), srtCol =45)
heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=colorRampPalette(c("#F7FBFF", "#08306B"))(n = 50),margins=c(12,12), srtCol =45 ,cexRow = 0.8, cexCol = 0.8)

#######################
#	NON-stringent DEPLETED
#######################

df_CpcogFull = cbind()
COG_distribution_all = data.frame(cbind(Var1=as.character(unique(ANVIO$COG_CATEGORY_ACC))))
for(group in unique(ANVIO_cat$group)){
	scoary <- read.table(paste('~/v4_SCOARY/',group,'_02_03_2018_1222.results.csv',sep=''),header=TRUE,sep=",") #protein clusters
	scoary = right_join(DEDUPLICATED3, scoary,by=('genes'))
	scoary = scoary[order(scoary$Benjamini_H_p),]
	scoary = scoary[scoary$Naive_p <= 0.01,]
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
}

names(df_CpcogFull)= c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
rownames(df_CpcogFull)=as.character(unique(ANVIO_cat$group))
df_CpcogFull = data.frame(t(df_CpcogFull))
df_CpcogFull$Var1 = c( "C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )

df_CpcogFull = right_join(COG_meaning, df_CpcogFull,by=('Var1'))
df_CpcogFull2 = df_CpcogFull[,3:ncol(df_CpcogFull)]
rownames(df_CpcogFull2) = paste(df_CpcogFull$Var1 , df_CpcogFull$long, sep=' - ')

heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=viridis(50),margins=c(12,12), srtCol =45)
heatmap.2(as.matrix(df_CpcogFull2),trace="none",col=colorRampPalette(c("#FFF5F0", "#67000D"))(n = 50),margins=c(12,12), srtCol =45,cexRow = 0.8, cexCol = 0.8)



#######################

DEDUPLICATED3 =  DEDUPLICATED[,c("protein_cluster_id",'COG_CATEGORY_ACC','COG_CATEGORY','COG_FUNCTION_ACC','COG_FUNCTION')]
DEDUPLICATED3$COG_FUNCTION = sub('\\\t.*', '', DEDUPLICATED3$COG_FUNCTION)
names(DEDUPLICATED3)= c("genes","COG_CATEGORY_ACC","COG_CATEGORY","COG_FUNCTION_ACC",'COG_FUNCTION')

for(group in unique(ANVIO_cat$group)){
	scoary <- read.table(paste('~/v4_SCOARY/',group,'_02_03_2018_1222.results.csv',sep=''),header=TRUE,sep=",") #protein clusters
	scoary = right_join(DEDUPLICATED3, scoary,by=('genes'))
	write.csv(scoary , file= paste('~/v4_SCOARY/FUNCTION_',group,'_added.results.csv',sep=''),row.names = FALSE)
}


#######################

write.csv(Tab , file= paste('~/Tab','_added.results.csv',sep=''),row.names = FALSE)
write.csv(Mat , file= paste('~/Tab','_added.results.csv',sep=''),row.names = FALSE)



names(COG_distribution_all)=c('Var1',as.character(unique(ANVIO_cat$group)))
#COG_distribution_all
COG_distribution_all = right_join(COG_meaning, COG_distribution_all,by=('Var1'))
COG_distribution_all$combined = paste(COG_distribution_all $Var1, COG_distribution_all $long)
COG_distribution_all2 =  COG_distribution_all[,3:11]
rownames(COG_distribution_all2)= COG_distribution_all $combined
COG_distribution_all2 = COG_distribution_all2[rowSums(is.na(COG_distribution_all2))!=ncol(COG_distribution_all2), ]
COG_distribution_all2[is.na(COG_distribution_all2)] <- 0


distance.col = dist(t(as.matrix(COG_distribution_all2)), method = "euclidean")
cluster.col = hclust(distance.col, method = "average")
distance.row = dist(as.matrix(COG_distribution_all2), method = "euclidean")
cluster.row = hclust(distance.row, method = "average")

COG_distribution_all2[COG_distribution_all2 == 0] <- NaN

heat2 = heatmap.2( as.matrix(COG_distribution_all2),
           col = inferno(75),Rowv= as.dendrogram(cluster.row),srtCol=45,
           trace = "none", 
           na.color="grey90",
           cexRow = 0.5, cexCol = 0.6,Colv=as.dendrogram(cluster.col),
           margins=c(12,12),keysize=0.75)

heat2 = heatmap.2( as.matrix(t(COG_distribution_all2)),
           col = inferno(75),Rowv= as.dendrogram(cluster.col),srtCol=45,
           trace = "none", 
           na.color="grey90",
           cexRow = 0.5, cexCol = 0.6,Colv=as.dendrogram(cluster.row),
           margins=c(12,12),keysize=0.75)






COG_list = c("D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )
sum=0
mat = colSums(t(table(droplevels(scoary[grep("C", scoary$COG_CATEGORY.y),c('COG_CATEGORY_ACC.x','COG_FUNCTION')]))))
CpcogFull = data.frame(sum(mat))
for(i in COG_list){
	mat = colSums(t(table(droplevels(scoary[grep(i, scoary$COG_CATEGORY.y),c('COG_CATEGORY_ACC.x','COG_FUNCTION')]))))
	CpcogFull = cbind(CpcogFull, sum(mat))
}
names(CpcogFull)= c( "C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X" )










DEDUPLICATED3 =  DEDUPLICATED[,c("protein_cluster_id",'COG_CATEGORY_ACC','COG_CATEGORY','COG_FUNCTION_ACC','COG_FUNCTION')]
DEDUPLICATED3$COG_FUNCTION = sub('\\\t.*', '', DEDUPLICATED3$COG_FUNCTION)

names(DEDUPLICATED3)= c("genes","COG_CATEGORY_ACC","COG_CATEGORY","COG_FUNCTION_ACC",'COG_FUNCTION')
DF5 = right_join(DEDUPLICATED3, testreverse,by=('genes'))
DF5$COG_FUNCTION = sub('\\\t.*', '', DF5$COG_FUNCTION)


testreverse = cbind("genes"= c('PC_00003891','PC_00005937','PC_00005738','PC_00010398','PC_00004477','PC_00008049','PC_00005445','PC_00004613','PC_00003641','PC_00003231','PC_00003985','PC_00003974','PC_00003749','PC_00003625','PC_00002465','PC_00002555','PC_00005315','PC_00002857','PC_00005578','PC_00003951','PC_00006127','PC_00000602','PC_00001315','PC_00003828','PC_00003411','PC_00003319','PC_00005726','PC_00008517','PC_00007420','PC_00002770','PC_00002806','PC_00004083','PC_00007389','PC_00004057','PC_00003499','PC_00001104','PC_00000355','PC_00003162','PC_00002989','PC_00002925','PC_00002728','PC_00007355','PC_00008764','PC_00005408','PC_00000509','PC_00003569','PC_00004296','PC_00003453','PC_00000725'))




######Soil versus Root######
rm(Map,Tab,trait,df_traits,df_Tab,df_Tab_vals,tree,root_node,subtree)
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Map<-read.table("metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Tab<-read.table(mat,header=T,sep="\t",row.names=1,check.names=F)
Tab<-t(Tab)
#The logic is to subset first using the columns of the matrix
Map<-Map[match(rownames(Tab),Map$taxon_oid),]
#Subset soil and root at the same time
Map_soil<-as.character(Map$taxon_oid[which(Map$Classification=="soil")])
Map_ra<-as.character(Map$taxon_oid[which(Map$Root.NotRoot=="RA")])
myids<-c(Map_soil,Map_ra)
Map<-Map[match(myids,Map$taxon_oid),]
Map<-droplevels(Map)
Tab<-Tab[match(Map$taxon_oid,rownames(Tab)),]
#Remove cols that sum to zero
Tab<-Tab[,which(colSums(Tab)!=0)]
Tab<-t(Tab)

#Detine the trait RA in this case should be 1 and soil 0
trait<-c(rep(0,length(Map_soil)),rep(1,length(Map_ra)))
df_traits<-data.frame(TraitY=trait)
rownames(df_traits)<-Map$taxon_oid

df_Tab<-data.frame(genes=rownames(Tab))
df_Tab_vals<-as.data.frame(Tab)
df_Tab<-cbind(df_Tab,df_Tab_vals)

#Read the general phylogenetic tree
#tree<-read.tree("/pine/scr/i/s/isai/3837_march_2017/3837_genomes_31scg_june2016.newick")
tree<-read.tree("3837_genomes_31scg_june2016.newick")
#Root the tree first using MRCA
#Use the firmicutes
#Paenibacillus Isolate 2517572151
#Bacillus Isolate 2623620997
root_node<-phytools::findMRCA(tree=tree,tips=c("2517572151","2623620997"))
tree<-reroot(tree,node.number=root_node)

subtree<-drop.tip(phy=tree,tip=which(!(tree$tip.label%in%Map$taxon_oid)))


write.table(df_traits,"traits_scoary_ra_soil.tsv",row.names=T,col.names=NA,sep=",",append=F,quote=F)
write.table(df_Tab,"matrix_scoary_ra_soil.tsv",row.names=F,col.names=T,sep=",",append=F,quote=F)
write.tree(subtree,"tree_scoary_ra_soil.newick")




#PlotHeatmap_res <- function(tree,Tab,Map,group){
 	#COGselect = head(finRes2$COG_FUNCTION_ACC,100)
    #COGselect = as.character(COGselect)
    #COGselect = COGselect[COGselect!='']
    #mat = t(table(droplevels(subset(ANVIO, COG_FUNCTION_ACC  %in% COGselect)[,c('COG_FUNCTION','genome_name')])))

	#dmat = data.frame(mat)
	#dmat = reshape(dmat, idvar = "genome_name", timevar = "COG_FUNCTION", direction = "wide")

	#colnames(dmat) = gsub("Freq.",'',colnames(dmat))
	#rownames(dmat) = dmat $genome_name
	#dmat = dmat[,2:length(colnames(dmat))]
	#dmat$genome = rownames(dmat)
	#matFull = right_join(dmat,dummy, by=c('genome'))
	#matFull[is.na(matFull)]=0
	#rownames(matFull) = matFull$genome
	#matFull = matFull[,1:length(colnames(matFull))-1]
	#mat = matFull

    #mat_HEAT = mat[,1:ncol(mat)]
	#mat_HEAT2 = mat_HEAT
	##mat_HEAT2[mat_HEAT2>1] <-1
	#mat_HEAT2_ordered <- mat_HEAT[new_order,]

	#distance.col = dist(t(as.matrix(mat_HEAT2_ordered)), method = "euclidean")
	#cluster.col = hclust(distance.col, method = "average")
	#mat_HEAT2_ordered[mat_HEAT2_ordered == 0] <- NaN

	#pdf(paste(group,"enriched_vs_others_heat.pdf",sep=''))
	#heat2 = heatmap.2( as.matrix(mat_HEAT2_ordered),
    #       col = inferno(75),Rowv= rep_tree_d,srtCol=45,
    #       trace = "none", 
    #       na.color="grey90",
    #       cexRow = 0.5, cexCol = 0.6,Colv=as.dendrogram(cluster.col),
    #       margins=c(12,12),keysize=0.75,main= paste("Enriched ",group, " versus other"))
	#dev.off()    
#}



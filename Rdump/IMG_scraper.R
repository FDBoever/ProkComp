#
# 	Frederik De Boever, 10-Jun-2018
#
# 	IMG_Scraper.R
#	Scraping the surface of IMG abundance tables using IMG_Scraper.R, can use PFAM, TIGRFAM, COG, or any other IMG category
#	It is designed to compare differences within genera as well as between genera
#
#	Recommended appraoch; 
#	1) 	go to img.jgi.doe.gov; genome search; and search for keywords; example; Oceanospirilalles
#	2) 	select all; and add to genome cart; repeat if you want to include more genomes: (you can select many different genera)
#	3) 	go to genome cart; scroll to bottom and select more fields to include in the table; select all genomes and export the taxon table
#	4)	For abundance tables; go to "Compare Genomes", "Abundance Profiles", "Function Comparison", select the genomes of interest
#	5)	run the analysis on any of the functions; for example, COG, and wait untill it's computed
#	6) 	select all, and click "Download tab-delimited file for Excel" 
#
#	This Script is designed to 
#	1) 	visual inspection of IMG abundance tables
#	2) 	combine meta data and abundance data to explore patterns
#	3) 	exploring subsets of genes using grep
#	4)	Quickly visualise ordination diagrams (genomes in functional space)
#	5)	calculate correlation between genomes on the basis of functional profiles
#
#	Reference:
#	IMG website: img.jgi.doe.gov


#########################################
#			LOAD packages
#########################################
library(gplots)
library(ggbiplot)
library(ggtree)
library(vegan)
library(reshape)
library(ape)
print('::: loading packages ::: done')


#########################################
#			INPUT
#########################################
# DOWNLOAD THE ABUNDANCE TABLES FROM IMG, AND IMPORT THEM HERE (needs raw abundance, not COG category)
#blaa = read.table("~/Desktop/documenten/testGamma.txt",header=TRUE)
#abundance_img = read.csv("~/Downloads/abundComp_pfam_27257.tab.xls",header=TRUE,sep="\t")
#abundance_img = read.csv("~/Downloads/abundComp_cog_5751.tab.xls",header=TRUE,sep="\t")
#abundance_img = read.table("~/Desktop/documenten/testGamma.txt",header=TRUE)
abundance_img = read.csv("~/Downloads/abundComp_cog_118403.tab.xls",header=TRUE,sep="\t")

# DOWNLOAD THE TAXON TABLE FROM IMG, AND IMPORT THEM HERE (will serve as metadata, remains optional)
metadata_img = read.csv("~/Downloads/taxontable118143_09-jun-2018.xls",header=TRUE,sep="\t")


# If you included a larger diversity in your abundance table, define your genus of interest
GenusOfInterest = "Marinobacter"

#Note Since is used to grep another '_' is added to prevent inclusion of similarly named genera, such as for example Marinobacter:Marinobacterium
GenusOfInterest = paste(GenusOfInterest, "_",sep='')

outdir = "~/IMG_SCRAPER/test/"
dir.create(outdir)
print('::: importing data ::: done')
print('::: making directory ... ')
print(paste('::: files will be stored in :::',outdir))

#########################################
#		SILVA BASED PHYLOGENY
#########################################
SilvaTree = read.tree('~/Documents/LTPs128_SSU_tree.newick')
selected = SilvaTree$tip.label[grepl(GenusOfInterest,SilvaTree$tip.label)]
subtree = drop.tip(SilvaTree,SilvaTree$tip.label[!grepl(GenusOfInterest,SilvaTree$tip.label)])
print('::: Loading Silva Tree ::: done')

## one could test this if needed

pdf(paste(outdir,paste('Silva_Genus_Of_Interest_unrooted',".pdf", sep = "_"), sep = "/"))
ggtree(subtree,layout = 'unrooted')+geom_treescale()
dev.off()
pdf(paste(outdir,paste('Silva_Genus_Of_Interest_rectangular',".pdf", sep = "_"), sep = "/"))
ggtree(subtree,layout = 'rectangular',ladderize=TRUE)+geom_treescale()
dev.off()
#pdf(paste(outdir,paste('Silva_Full_tree',i,".pdf", sep = "_"), sep = "/"))
#ggtree(SilvaTree,layout = 'unrooted')+geom_treescale()
#dev.off()





	
	

#########################################
#		IMG DATAFRAME CONSTRUCTION
#########################################
# there is no need to use excell anylonger, just run the below, and get the data ready to be mined
abundance_img  = abundance_img[,c(1,2,seq(4, length(colnames(abundance_img)), 3))]

tidyUp = c("_pvalue","_geneCount","_contamination_screened_")
for(i in tidyUp){colnames(abundance_img) = gsub(i, "", colnames(abundance_img))}

print('::: Data organisation ::: done')

#########################################
#		OPTIONAL STRAIN SELECTION
#########################################

#StrainList = c("algicola","VT8","hydrocarbonoclastiucs","EN3","MCT","EN3","SM19","DSM")
#toMatch = c("Func_name","Func_id", StrainList)
#abundance_img  = abundance_img[, grep(paste(toMatch,collapse="|"), colnames(abundance_img))]
#print('::: Strain selection ::: done')


############
# CAUTION; INSPECT THE GRAPHICS BELOW, TO FILTER OUT PUTITQTIVELY LOW QUQLITY GENOMES OUT

fullheatdf = abundance_img
rownames(fullheatdf)= fullheatdf$Func_id
fullheatdf= fullheatdf[,-c(1,2)]
Totals = data.frame(as.matrix(colSums(fullheatdf)),colnames(fullheatdf))
colnames(Totals)=c('predicted','genome')


Totals$genus = gsub("\\_.*","",Totals$genome)

maxim = c()
for(i in Totals$genus){maxim = c(maxim,max(Totals[Totals$genus == i,"predicted"]))}

Totals$maxim= maxim
hist(Totals$predicted)

Totals$putative= Totals$predicted/Totals$maxim
Totals$putative_threshold = ifelse(Totals$putative > 0.4, "good", "bad")
Totals[Totals$predicted < 1500,]$putative_threshold = "bad"

hist(Totals$predicted)

pdf(paste(outdir,paste('Predicted_genes_Filtering_bar_per_genome',i,".pdf", sep = "_"), sep = "/"))
ggplot(Totals, aes(x = genome,y= predicted,fill= putative_threshold)) + geom_bar(stat="identity")+coord_flip()+ theme_light()+theme(axis.line = element_line(size=0.5),axis.ticks = element_line(size=0.5,colour="black"),axis.text.x = element_text(colour="black",angle = 90, hjust = 1, size=7),axis.text.y = element_text(colour="black",size=7))
dev.off()

pdf(paste(outdir,paste('Predicted_genes_Filtering_points_per_genus',i,".pdf", sep = "_"), sep = "/"))
ggplot(Totals, aes(y= predicted,x=genus,color= putative_threshold)) + theme_light()+theme(axis.line = element_line(size=0.5),axis.ticks = element_line(size=0.5,colour="black"),axis.text.x = element_text(colour="black",angle = 90, hjust = 1, size=7),axis.text.y = element_text(colour="black",size=7))+geom_jitter(size=0.3)
dev.off()

pdf(paste(outdir,paste('Predicted_genes_Filtering_boxplot_per_genus',i,".pdf", sep = "_"), sep = "/"))
ggplot(Totals, aes(y= predicted,x=genus,color= putative_threshold)) + geom_boxplot()+geom_jitter(size=0.3) + theme_light()+theme(axis.line = element_line(size=0.5),axis.ticks = element_line(size=0.5,colour="black"),axis.text.x = element_text(colour="black",angle = 90, hjust = 1, size=7),axis.text.y = element_text(colour="black",size=7))
dev.off()


Totals[Totals$putative_threshold == "bad","genome"]


pdf(paste(outdir,paste('metadata_quality_vs_size_all',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(x= reorder(High.Quality, Genome.Size.....assembled),y= Genome.Size.....assembled),fill=Totals$putative_threshold) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))
dev.off()

ggplot(metadata_img, aes(x= reorder(Genus, Genome.Size.....assembled),y= Genome.Size.....assembled,fill=Totals$putative_threshold)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))

ggplot(metadata_img, aes(x= reorder(Genus, Genome.Size.....assembled),y= Genome.Size.....assembled,fill=High.Quality)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))


ggplot(metadata_img, aes(x= reorder(High.Quality, Genome.Size.....assembled),y= Genome.Size.....assembled)) + geom_boxplot()+geom_jitter(size=0.5,colour=Totals$putative_threshold	*1)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))


#########################################
#		FILTER IF NEEDED
#########################################
fullheatdf = fullheatdf[,c(as.character(Totals[Totals$putative_threshold == "good","genome"]))]
abundance_img = abundance_img[,c("Func_id","Func_name",as.character(Totals[Totals$putative_threshold == "good","genome"]))]
Totals = Totals[c(as.character(Totals[Totals$putative_threshold == "good","genome"])),]
metadata_img = metadata_img[metadata_img$Genome.Name...Sample.Name %in% gsub("_"," ",c(as.character(Totals[Totals$putative_threshold == "good","genome"]))),]


ggplot(Totals, aes(y= predicted,x=reorder(genus,predicted))) + geom_boxplot()theme_light()+theme(axis.line = element_line(size=0.5),axis.ticks = element_line(size=0.5,colour="black"),axis.text.x = element_text(colour="black",angle = 90, hjust = 1, size=7),axis.text.y = element_text(colour="black",size=7))

#		TEST TREES
#########################################

subtree = drop.tip(SilvaTree,SilvaTree$tip.label[!grepl(paste(Totals$genus,collapse="|"),SilvaTree$tip.label)])
print('::: Loading Silva Tree ::: done')

pdf(paste(outdir,paste('Silva_all_genera_in_data_unrooted',".pdf", sep = "_"), sep = "/"))
ggtree(subtree,layout = 'unrooted')+geom_treescale()
dev.off()
pdf(paste(outdir,paste('Silva_all_genera_in_data_rectangular',".pdf", sep = "_"), sep = "/"))
ggtree(subtree,layout = 'rectangular',ladderize=TRUE)+geom_treescale()
dev.off()

genera = unique(Totals$genus)
l <- list();	
	for(genus in genera){
		l[[genus]] = as.character(subtree$tip.label[grep(genus,subtree$tip.label,ignore.case=TRUE)])
	}
annotatedTree <- groupOTU(subtree, l)

#### PLOT TREES
pdf(paste(outdir,paste('SILVA_gzoom',".pdf", sep = "_"), sep = "/"))
gzoom(annotatedTree, grep(GenusOfInterest, annotatedTree$tip.label))
dev.off()
pdf(paste(outdir,paste('SILVA_circular_annotated_circ',".pdf", sep = "_"), sep = "/"))
ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,) +geom_treescale()+geom_tiplab(aes(color= group,angle=angle),size=2,align=TRUE,linesize=.2)
dev.off()
pdf(paste(outdir,paste('SILVA_circular_annotated_rect',".pdf", sep = "_"), sep = "/"))
ggtree(annotatedTree)+geom_tippoint(aes(color= group), alpha=1) +geom_treescale()+geom_tiplab(aes(color= group),size=2)
dev.off()
pdf(paste(outdir,paste('SILVA_circular_annotated_circ_nolab',".pdf", sep = "_"), sep = "/"))
ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,)
dev.off()


#### HIGHLIGHT CLADES
genusTip = annotatedTree$tip.label[grepl(paste(Totals$genus,collapse="|"),annotatedTree$tip.label)]
rca = getMRCA(annotatedTree,tip= genusTip)

p3.1 = ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,)
colint = 200
for(genus in unique(genera)){
	genusTip = annotatedTree$tip.label[grepl(paste(genus,'_',sep=''),annotatedTree$tip.label)]
	rca = getMRCA(annotatedTree,tip= genusTip)
	print(rca)
	p3.1 = p3.1 + geom_hilight(node= rca, fill=colors()[colint])+geom_cladelabel(node= rca, label= genus,offset.text=0.05)
	colint = colint + 10 
	}
pdf(paste(outdir,paste('SILVA_highlighted_clades',".pdf", sep = "_"), sep = "/"))
ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,)
dev.off()

#########################################
#		MONOPHYLY+
#########################################

monophyl = c()
for(genus in unique(genera)){
	monophyl = c(monophyl,is.monophyletic(phy = annotatedTree, tips = annotatedTree$tip.label[grepl(genus, annotatedTree$tip.label)]))
	}
monophyl = data.frame( cbind("genus" = genera, "monophyl" =  monophyl))

print(monophyl)

##########	annotate only monophyletic groups		###########

p3.2 = ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,)
colint = 200
for(genus in unique(genera)){
	if(is.monophyletic(phy = annotatedTree, tips = annotatedTree$tip.label[grepl(genus, annotatedTree$tip.label)])){
		genusTip = annotatedTree$tip.label[grepl(paste(genus,'_',sep=''),annotatedTree$tip.label)]
		rca = getMRCA(annotatedTree,tip= genusTip)
		print(rca)
		p3.1 = p3.1 + geom_hilight(node= rca, fill=colors()[colint])+geom_cladelabel(node= rca, label= genus,offset.text=0.05)
		colint = colint + 10 
	}
}
pdf(paste(outdir,paste('SILVA_only_monophyletic_clades_highlighted',".pdf", sep = "_"), sep = "/"))
print(p3.2)
dev.off()



#########################################
#		METADATA_SCANNER
#########################################
pdf(paste(outdir,paste('metadata_size_vs_culttype',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(y= Genome.Size.....assembled ,x= Culture.Type)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))
dev.off()

pdf(paste(outdir,paste('metadata_size_per_genus',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(y= Genome.Size.....assembled ,x= reorder(Genus, Genome.Size.....assembled))) + geom_boxplot()+geom_jitter(size=0.5)+theme_light()+theme(axis.line = element_line(size=0.5),axis.ticks = element_line(size=0.5,colour="black"),axis.text.x = element_text(colour="black",angle = 90, hjust = 1, size=7),axis.text.y = element_text(colour="black",size=7))+xlab("genus")+ylab("genome size (bp)")
dev.off()

pdf(paste(outdir,paste('metadata_size_per_family',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(y= Genome.Size.....assembled ,x= Family)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))
dev.off()

pdf(paste(outdir,paste('metadata_size_vs_GOLD.Sequencing.Strategy',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(y= Genome.Size.....assembled ,x= GOLD.Sequencing.Strategy)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))
dev.off()

pdf(paste(outdir,paste('metadata_size_vs_Habitat',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(y= Genome.Size.....assembled ,x= Habitat)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))
dev.off()

pdf(paste(outdir,paste('metadata_size_vs_COG.Count',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(x= Genome.Size.....assembled ,y= COG.Count.....assembled)) +geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))
dev.off()

pdf(paste(outdir,paste('metadata_size_vs_GC',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(x= Genome.Size.....assembled ,y= GC.....assembled)) +geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))
dev.off()


pdf(paste(outdir,paste('metadata_GC_vs_genus',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(x= reorder(Genus, GC.....assembled),y= GC.....assembled)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))
dev.off()

pdf(paste(outdir,paste('metadata_GC_vs_genus2',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(x= reorder(Genus, GC.....assembled),y= GC.....assembled)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))
dev.off()


pdf(paste(outdir,paste('metadata_quality_vs_size_filtered',".pdf", sep = "_"), sep = "/"))
ggplot(metadata_img, aes(x= reorder(High.Quality, Genome.Size.....assembled),y= Genome.Size.....assembled)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))
dev.off()

#########################################
#		FULL HEATMAP
#########################################

allcog = dim(fullheatdf)
redheatdf = fullheatdf[rowSums(fullheatdf)>0,]
redcog = dim(redheatdf)
print(paste("your data holds" , round((redcog/allcog)[1], digits=4) , "% of the total COG:PFam acessions in IMG"))

distance_col = dist(as.matrix(redheatdf),method="euclidean")
cluster_col = hclust(distance_col,method="average")
distance_now = dist(as.matrix(t(redheatdf)),method="euclidean")
cluster_now = hclust(distance_now,method="average")

redheatdf[redheatdf ==0] = NaN

pdf(paste(outdir,paste('FULL_HEATMAP',".pdf", sep = "_"), sep = "/"))
heatmap.2(as.matrix(t(redheatdf)),col = colorpanel(100,"black","yellow","green"),trace="none",na.color="grey90",cexRow=0.2, cexCol=0.2,margins=c(5,5),Rowv 	=as.dendrogram(cluster_now),Colv=as.dendrogram(cluster_col))
dev.off()

#########################################
#		SELECT A FUNCTION AND MINE THE DATASET
#########################################

#
FuncOfInterest = "transposase"

mat = abundance_img[grep(FuncOfInterest, abundance_img$Func_name),]
mat2 = mat[,3:length(names(mat))]
rownames(mat2) = paste(mat$Func_name,mat$Func_id,sep="-")
mat2 = mat2[rowSums(mat2)>0,]
DNF <- ifelse(colSums(mat2) > 1, "present", "not-present")
DNF  = data.frame(DNF)
DNF$genome  = rownames(DNF)

mat3 =  melt(t(mat2))
colnames(mat3) = c("genome","func","value")
combined = merge(x = mat3, y = DNF, by = "genome", all = TRUE)
ggplot(combined, aes(x = func, y = value)) + geom_boxplot()+coord_flip()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))+theme_bw()

Sumtotals = data.frame(as.matrix(colSums(mat2)),colnames(mat2))
colnames(Sumtotals)=c('FuncSelect','genome')

Sumtotals$genus = gsub("\\_.*","", Sumtotals$genome)
Sumtotals$DNF = DNF $DNF
ggplot(Sumtotals, aes(x = genome,y= FuncSelect)) + geom_bar(stat="identity")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))

ggplot(Sumtotals, aes(x = genus,y= FuncSelect)) + geom_bar(stat="identity")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))

ggplot(Sumtotals, aes(x = genus,y= FuncSelect)) +geom_point()+coord_flip()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))+geom_jitter()

ggplot(Sumtotals, aes(x = genus,y= DNF)) +geom_point()+coord_flip()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))+geom_jitter()

ggplot(Sumtotals, aes(x = genus,y= DNF)) + geom_bar(stat="identity")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))

ggplot(Sumtotals, aes(x = reorder(genus, FuncSelect),y= FuncSelect,fill= DNF)) + geom_bar(stat="identity")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))

#########################################
#		SELECTIVE HEATMAPS
#########################################

SelectHeat = function(toMatch){
	mat = abundance_img[grep(paste(toMatch,collapse="|"), abundance_img$Func_name,ignore.case=TRUE),]
	mat2 = mat[,3:length(names(mat))]
	rownames(mat2) = paste(mat$Func_name,mat$Func_id,sep="-")
	mat2=data.frame(mat2)
	mat2 = mat2[rowSums(mat2)>0,]
	#mat2 = mat2[-which(!unlist(lapply(mat2,function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))),]

	#mat3 =  melt(t(mat2))
	#colnames(mat3) = c("genome","func","value")
	#combined = merge(x = mat3, y = DNF, by = "genome", all = TRUE)
	#ggplot(combined, aes(x = func, y = value)) + geom_boxplot(aes(fill=DNF))+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6))

	distance_col = dist(as.matrix(mat2),method="euclidean")
	cluster_col = hclust(distance_col,method="average")
	distance_now = dist(as.matrix(t(mat2)),method="euclidean")
	cluster_now = hclust(distance_now,method="average")

	mat2[mat2==0] = NaN
	
	heatmap.2(as.matrix(t(mat2)),col = colorpanel(100,"black","yellow","green"),trace="none",na.color="grey90",cexRow=0.2, cexCol=0.2, Rowv 	=as.dendrogram(cluster_now),Colv=as.dendrogram(cluster_col),margins=c(10,10),main=paste(toMatch,sep=" "))

}

#########################################
# Plot the heatmaps for selective functions; examples below

SelectHeat(c("permease"))
SelectHeat(c("secretion"))
SelectHeat(c("flagel"))
SelectHeat(c("transposase","integrase"))

#########################################


#make lists of functions you want to visualise
#THE BELOW EXAMPLE IS FOR COG BASED ABUNDANCE TABLE

GrepSearches = list(c("transport","nitr","nitrate","sensing","efflux","sulfate","iron","NAD","FAD","quinone"),c("NAD","FAD","quinone"),c("nitr","nitrate","sensing","efflux","sulfate","iron","NAD","FAD","quinone"),c("Bacterial","regulatory"),c("sodium","Sodium","Na+"),c("biotin","Biotin"),c("CoA","pyruvate","Pyruvate"),c("ose"),c("channel","Channel"),c("oxidoreductase"),c("oxidase"),c("transferase"),c("methionine"),c("carboxyl"),c("phosphatase"),c("hydrolase") ,c("redoxin") ,c("desaturase","fatty") ,c("synthase")  ,c("kinase")  ,c("domain")  ,c("protease")  ,c("flag")    ,c("redoxin")    ,c("toxin"),c("reductase") ,c("Dioxygenase","oxygenase")  ,c("binding")  ,c("gluta","Gluta")  ,c("transketo"),c("NAD","FAD"),c("quinone"))

for(i in GrepSearches){
	tryCatch({
	pdf(paste(outdir,paste('GrepBased',i,".pdf", sep = "_"), sep = "/"))
	SelectHeat(i)
	dev.off()
  	#Closure of tryCatch({ which prevents ERROR's to stop the loop
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#########################################
#		SIMPLE ORDINATION
#########################################

mat2 = fullheatdf[rowSums(fullheatdf)>0,]

mat2 = mat2[-which(!unlist(lapply(mat2,function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))),]
data.pca <- prcomp(t(mat2), scale. = TRUE)
g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, ellipse = FALSE, circle = FALSE,groups=Totals$genus, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')  


#########################################
#		Correlation graphic
#########################################

metaTax = metadata_img[,c("Genus","Family","Order")]   
metaTax = metaTax[!duplicated(metaTax),]
metaTax$genusX1 = metaTax$Genus
metaTax$genusX2 = metaTax$Genus

#if you need to reclassify
levels(metaTax$Family)=c(levels(metaTax$Family),"Marinobacteraceae")
metaTax[metaTax$Genus=="Marinobacter",c("Order")]= "Oceanospirillales"
metaTax[metaTax$Genus=="Marinobacter",c("Family")]= "Marinobacteraceae"

meltTax = function(corMat, metaTax){
	moltencorMat = melt(corMat)
	moltencorMat$species1 = gsub("^([^_]*_[^_]*)_.*$", "\\1", moltencorMat$X1)
	moltencorMat$species2 = gsub("^([^_]*_[^_]*)_.*$", "\\1", moltencorMat$X2)
	moltencorMat$genusX1 = gsub("\\_.*","", moltencorMat$X1)
	moltencorMat$genusX2 = gsub("\\_.*","", moltencorMat$X2)
	moltencorMat$genus_of_interest = ifelse(paste(moltencorMat$genusX1,"_",sep="")==GenusOfInterest,"yes","no")

	moltencorMat = merge((moltencorMat), metaTax[,c("genusX1", "Family", "Order")], by = 'genusX1')
	moltencorMat = merge((moltencorMat), metaTax[,c("genusX2", "Family", "Order")], by = 'genusX2')

	moltencorMat$taxorder = ifelse(moltencorMat$X1==moltencorMat$X2,"itself",ifelse(moltencorMat$species1==moltencorMat$species2 && moltencorMat$species1 !=paste(moltencorMat$genusX1 ,"sp",sep="_"),"intraspecies",ifelse(moltencorMat$genusX1==moltencorMat$genusX2,"intragenus",ifelse(moltencorMat$Family.x ==moltencorMat$Family.y,"intraFamily",ifelse(moltencorMat$Order.x == moltencorMat$Order.y,"intraOrder","intrerOrder")))))

	#ggplot(moltencorMat, aes(x= reorder(taxorder,value),y= value,colour= genus_of_interest)) +geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))+scale_colour_manual(values=c("black","orange"))

	return(moltencorMat)
}


redabundance = fullheatdf[rowSums(fullheatdf)>0,]

corMat = cor(redabundance,method="spearman")
heatmap.2(as.matrix(corMat),col = colorpanel(100,"darkblue","white","darkred"),trace="none",cexRow= 0.1, cexCol= 0.1,main="spearman")
data.moltencorMat =meltTax(corMat, metaTax)

corMat = cor(redabundance,method="pearson")
heatmap.2(as.matrix(corMat),col = colorpanel(100,"darkblue","white","darkred"),trace="none",cexRow= 0.1, cexCol= 0.1,main="pearson")
tmp.moltencorMat2 = meltTax(corMat, metaTax)
data.moltencorMat$Pearson = tmp.moltencorMat2$value

ggplot(data.moltencorMat, aes(x= reorder(taxorder,Pearson),y= value,colour= genus_of_interest)) +geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))+scale_colour_manual(values=c("black","orange"))

ggplot(data.moltencorMat, aes(x= Pearson,y= value,colour= genus_of_interest)) +geom_point(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))+scale_colour_manual(values=c("black","orange"))


ggplot(moltencorMat, aes(x= reorder(taxorder,value),y= value,colour= genus_of_interest)) +geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))+scale_colour_manual(values=c("black","orange"))



ggplot(moltencorMat, aes(x= reorder(taxorder,value),y= value,colour= genus_of_interest)) +geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))+scale_colour_manual(values=c("black","orange"))

ggplot(moltencorMat, aes(x= taxorder,y= value,colour= genus_of_interest)) +geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))+scale_colour_manual(values=c("black","orange"))



ggplot(moltencorMat, aes(x= taxorder,y= value,colour= genus_of_interest)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))




gsub("^([^-]*-[^-]*).*", "\\1", moltencorMat$X2)

ggplot(moltencorMat, aes(x= taxorder,y= value,fill)) + geom_boxplot()+geom_jitter(size=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=5))




           
#########################################
#		SIGNIFICANT FUNC
#########################################  
abund_tab  = t(fullheatdf)
groups = as.factor(Totals$genus)           

#Apply normalisation
data = abund_tab/rowSums(abund_tab)
data = data.frame(data)
kruskal.wallis.alpha = 0.01

kruskal.wallis.table <- data.frame()
for (i in 1:dim(data)[2]) {
  ks.test <- kruskal.test(data[,i], g=groups)
  # Store the result in the data frame
  kruskal.wallis.table <- rbind(kruskal.wallis.table,
                                data.frame(id=names(data)[i],
                                           p.value=ks.test$p.value
                                ))
  # Report number of values tested
  cat(paste("Kruskal-Wallis test for ",names(data)[i]," ", i, "/", 
            dim(data)[2], "; p-value=", ks.test$p.value,"\n", sep=""))
}

kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
 
kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, 
                                    size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
 
kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value,
                                                   decreasing=FALSE), ]
kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor

last.significant.element <- max(which(kruskal.wallis.table$q.value <= kruskal.wallis.alpha))
selected <- 1:last.significant.element
diff.cat.factor <- kruskal.wallis.table$id[selected]
diff.cat <- as.vector(diff.cat.factor)

df<-NULL
for(i in diff.cat){
  tmp<-data.frame(data[,i],groups,rep(paste(i," q = ",round(kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"],5),sep=""),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Type","Taxa")
 
p<-ggplot(df,aes(Type,Value,colour=Type))+ylab("Log-relative normalised")
p<-p+geom_boxplot()+geom_jitter()+theme_bw()+
  facet_wrap( ~ Taxa , scales="free", ncol=3)
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
print(p)
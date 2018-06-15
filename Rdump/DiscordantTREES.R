#DiscordantTREES.R


library(kdetrees)
library(ape)
library("seqRFLP")
library(msa)


#This is currently tailored to use ANVIO data and derive the Single Copy Orthologs (SCO) from there
GenomestoInclude = CheckMANVIO2$Bin_Id

outdir = "~/testFasta_2"
dir.create(outdir)

#remove all the genes unique to a single genome
COREDF = PA_AnvioData[rowSums(PA_AnvioData)== length(colnames(PA_AnvioData)), ]
COREDF = AnvioData[rownames(COREDF),]
COREDF = COREDF[rowSums(COREDF)== length(colnames(COREDF)),]
CORE_PCs = rownames(COREDF)

nr_SCO = dim(COREDF)[1]
COREanv = ANVIO[ANVIO$protein_cluster_id %in% CORE_PCs, ]
COREanv = COREanv[COREanv $genome_name %in% GenomestoInclude, ]


#-------------------------------------------------------------------------------
# 	KDE TREES
#-------------------------------------------------------------------------------
#	
treelist = list()
class(treelist) <- "multiPhylo"
i = 0
for(PC in CORE_PCs){
	i=i+1
	selectedCOREanv = COREanv[COREanv$protein_cluster_id == PC,]
	fasta_frame  = selectedCOREanv[,c('genome_name','aa_sequence')]
	df.fasta = dataframe2fas(fasta_frame,file=paste(outdir,paste(PC,".fasta", sep = ""), sep = "/"))

	sequences <- Biostrings::readAAStringSet(paste(outdir,paste(PC,".fasta", sep = ""), sep = "/"),format='fasta')

	aln = read.alignment(file = paste(outdir,paste(PC,".fasta", sep = ""), sep = "/"), format = "fasta")
	
	# dirty but fast: distance matrix and NJ tree reconstruction, and assesment of monophyly of a genus (CHANGE MARINOBACTER); 
	dist.aln <- dist.alignment(aln, "identity")
	njTree <- njs(dist.aln )
	njTreeMR = ladderize(midpoint.root(njTree), right = FALSE)
	print(paste(PC , "::", i,"/",length(CORE_PCs),"::",'tree building done',sep=''))
	
	treelist[[PC]]= njTreeMR
}
    
kdeRes =  kdetrees(treelist)
kdeRes2 = kdetrees(treelist, k=1.25, distance="geodesic", topo.only=FALSE)
kdeRes3 = kdetrees(treelist, k=1.25, distance="diss", topo.only=TRUE)
kdeRes4 = kdetrees(treelist, k=1.25, distance="diss", topo.only=FALSE)
kdeRes5 = kdetrees(treelist, k=1.25, distance="geodesic", topo.only=TRUE)

plot(kdeRes5)
hist(kdeRes5)
 
 
write.tree(treelist, paste(outdir,paste('multiPhylo',".tre", sep = ""), sep = "/"), tree.names = TRUE)
write.tree(kdeRes5$outliers,file=paste(outdir,paste('outliers',".tre", sep = ""), sep = "/"))
result.df <- as.data.frame(kdeRes5)
write.csv(result.df, file=paste(outdir,paste('KDE_Scores',".csv", sep = ""), sep = "/"))


outlierGenes = COREanv[COREanv$protein_cluster_id %in% names(kdeRes5$outliers), ]
coalescentGenes = COREanv[!(COREanv$protein_cluster_id %in% names(kdeRes5$outliers)), ]

outlierGenesDEDUPLICATED =  subset(outlierGenes, !duplicated(protein_cluster_id))[,c('protein_cluster_id','COG_CATEGORY_ACC','COG_FUNCTION_ACC', 'COG_FUNCTION')]
coalescentGenesDEDUPLICATED =  subset(coalescentGenes, !duplicated(protein_cluster_id))[,c('protein_cluster_id','COG_CATEGORY_ACC','COG_FUNCTION_ACC', 'COG_FUNCTION')]

write.csv(outlierGenesDEDUPLICATED, file=paste(outdir,paste('outlier_SCO',".csv", sep = ""), sep = "/"))
write.csv(coalescentGenesDEDUPLICATED, file=paste(outdir,paste('coalescent_SCO',".csv", sep = ""), sep = "/"))


#	COPY THE FASTA-FILES TO A DISTINCT FOLDER
#-------------------------------------------------------------------------------

coalescentDIR = paste(outdir,'/COALESCENT',sep='')
dir.create(coalescentDIR)

fasta_list <- list.files(outdir, "\\.fasta$",full.names=T)
fasta_list = fasta_list[grepl(paste(unique(as.character(coalescentGenes$protein_cluster_id)),collapse='|'), fasta_list)]

# copy the files to the new folder
file.copy(fasta_list, coalescentDIR)



#-------------------------------------------------------------------------------
# 	DISTANCE BASED
#-------------------------------------------------------------------------------
#	

ANVIO_cat3 = ANVIO_cat2
rownames(ANVIO_cat3)= ANVIO_cat3$Bin_Id
ANVIO_cat3 = ANVIO_cat3[colnames(AccesoryDF),]



for(tree in kdeRes5$outliers){
	cophTree = cophenetic(tree)
	ph_cophTree = otu_table(cophTree,taxa_are_rows=TRUE)
	
	sampledata = sample_data(ANVIO_cat3[tree$tip.label,])
	sample_names(sampledata) = sample_names(ph_cophTree)
	phyloseq = phyloseq(ph_cophTree, sampledata)

	iMDS  <- ordinate(phyloseq, "NMDS", distance='jaccard')
	p <- plot_ordination(phyloseq, iMDS, color="group2")
	# Add title to each plot
	p <- p + geom_point(size=1) + scale_colour_manual(values = c("#39811D","#86DA81","#A6D8D4","#F2A968","#F2EC70","#E38FDD","#898989","#76AECF","#B34D22"))+ theme(legend.position="none",axis.line = element_line(colour="black"),axis.ticks = element_line(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())  + ggtitle(paste("ourlier: ", sep=""))
	print(p)
	
	hist(melt(cophTree)$value)

}
	



	

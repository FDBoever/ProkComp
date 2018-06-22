##############################################################################################################################################
#	Last Edited: 14-03-2018
#
#	BLASTp.R
#
# 	This script takes a fasta file containing single or multiple AA sequences and will for each sequence:
#		1) query NCBI using blastp
#		2) multi-sequence alignment of blasthits using CLUSTAL 2.1
#		3) infer a quick Neighbor Joining phylogenetic tree and visualise
#		4) check whether sequences retrieved from your genus of interest are forming a monophyletic clade
#		5) store taxonomic information of blast hits
#		
#	INPUT:
#		- InputFasta: 			in code below, set path to your fasta file
#		- GenusOfInterest		specify genus of interest (for monophyly detection and more)
# 		- outdir			 	fill in your prefered output directory (no "/" at the end!!!) 
#		- TaxTable		 		change to the path directing to the TaxTable file 
#								TaxTable is derived from recent classification page of Bacterio.net (late 2017)
#
#	OUTPUT per sequence:
#		- blast output in tab seperated .txt file
#		- fasta files, both unaligned as CLUSTAL aligned (consider to use these to build your own trees)
#		- pdfs with graphic representation of quick and dirty trees 
#
#
##############################################################################################################################################

#-------------------------------------------------------------------------------
# 	install and load packages
#-------------------------------------------------------------------------------
#Note, some packages need to be installed via Bioconductor 
#source("https://bioconductor.org/biocLite.R")
#biocLite("annotate")
#biocLite("Biostring")
#biocLite("msa")
#biocLite("tidyr")
#biocLite("bios2mds")
#biocLite("ggtree")

#install.packages(c("seqRFLP"),dependencies=TRUE)

library(seqinr)
library("seqRFLP")
library(msa)
library(ape)
library(annotate)
library(phytools)
library(stringr)
library(Biostring)
library(plyr)
library(bios2mds)
library(ggtree)

#-------------------------------------------------------------------------------
# 	IGNORE, this is how it is tied up to ProkComp
#-------------------------------------------------------------------------------

#genes_of_interest = c('PC_00005278','PC_00005173','PC_00005166','PC_00005288','PC_00005191','PC_00005167','PC_00005226','PC_00005295','PC_00005307','PC_00005287','PC_00005346','PC_00005242','PC_00005368','PC_00005190','PC_00005168','PC_00005201','PC_00005351','PC_00005397','PC_00005284','PC_00004976','PC_00005029','PC_00005080','PC_00005088','PC_00004947')
#genes_of_interest = selectedSums
#genes_of_interest = as.character(DEDUPLICATED[DEDUPLICATED $protein_cluster_id %in% selectedSums & DEDUPLICATED$COG_CATEGORY_ACC %in% c('C'), ]$protein_cluster_id)


#-------------------------------------------------------------------------------
#  INPUT
#-------------------------------------------------------------------------------
# - InputFasta			change to the path directing to a fasta file containing your sequences of interest
# - GenusOfInterest 	fill in your genus of interest, 
# - outdir			 	fill in your prefered output directory (no "/" at the end!!!) 
# - TaxTable		 	change to the path directing to the TaxTable file 


InputFasta = readAAStringSet('~/nar.fasta',format='fasta')
GenusOfInterest = "Marinobacter"
outdir = "~/nar2"
TaxTable = read.table("~/TaxTable.txt",sep='\t')

#Note Since is used to grep another space is added to prevent inclusion of similarly named genera, such as for example Marinobacter:Marinobacterium
GenusOfInterest = paste(GenusOfInterest, " ",sep='')

dir.create(outdir)

#SilvaTree = read.tree('LTPs128_SSU_tree.newick')


#-------------------------------------------------------------------------------
# 	MULTIPLOT FUNCTION TO MAKE COMPOSITE GRAPHICS
#-------------------------------------------------------------------------------

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


#-------------------------------------------------------------------------------
# 	MAIN CODE
#-------------------------------------------------------------------------------
# Runs over all provided sequences

FastaFiles <- list.files(path="~/Genomics/FASTA", pattern="*.fasta", full.names=T, recursive=FALSE)


#Note Since is used to grep another space is added to prevent inclusion of similarly named genera, such as for example Marinobacter:Marinobacterium
GenusOfInterest = paste(GenusOfInterest, " ",sep='')

dir.create(outdir)
~/Genomics/NifH.fasta
for(file in FastaFiles){
	InputFasta = readAAStringSet(file,format='fasta')
	GenusOfInterest = "Marinobacter "
	outdir = paste('~/',substring(basename(file),1,nchar(basename(file))-6),sep='')
	dir.create(outdir)
	print(paste('STARTING ANALYSIS FOR',basename(file)))
	blastpFASTA(InputFasta, GenusOfInterest, outdir)
}

InputFasta = readAAStringSet('~/Genomics/NifH.fasta',format='fasta')
GenusOfInterest = "Marinobacter "
outdir ='~/NifH'
dir.create(outdir)


#-------------------------------------------------------------------------------
# 	FROM FASTA
#-------------------------------------------------------------------------------


blastpFASTA <- function(InputFasta, GenusOfInterest, outdir ){
OverallTaxonData = rbind()
topHit = c()

#for(PC in genes_of_interest){
for(seqnr in 1:length(InputFasta)){
	PC = InputFasta@ranges@NAMES[seqnr]
	tryCatch({
	print(paste("Currently starting the analysis for:",PC,sep=' '))
	#res = blastSequences(as.character(ANVIO[ANVIO$protein_cluster_id== PC,'aa_sequence'][1]), "nr", program="blastp", 	as="XML",hitListSize="100",timeout=2000)
	
	#performing the blast search, generally slow, (tweak maxium waiting time by changing 'timeout' value )
	res = blastSequences(as.character(InputFasta[seqnr]), "nr", program="blastp", as="XML",hitListSize="3",timeout=2000)
	print('res done')
		
	# blastSequences returns output in XML format, lukily it is convenently structured 
	blast_out = data.frame(cbind(cbind(	
		Hit_num = sapply(res["//Hit_num"], xmlValue),
		Hit_id= sapply(res["//Hit_id"], xmlValue), 
		Hit_def = sapply(res["//Hit_def"], xmlValue),
		Hit_accession = sapply(res["//Hit_accession"], xmlValue),
		Hit_len= sapply(res["//Hit_len"], xmlValue)),
		cbind(	
		Hsp_bit_score= sapply(res["//Hsp_bit-score"], xmlValue), 
		Hsp_score = sapply(res["//Hsp_score"], xmlValue),
		Hsp_evalue = sapply(res["//Hsp_evalue"], xmlValue), 
		Hsp_identity = sapply(res["//Hsp_identity"], xmlValue),
		Hsp_gaps = sapply(res["//Hsp_gaps"], xmlValue), 
		Hit_len= sapply(res["//Hit_len"], xmlValue),
		Hsp_align_len = sapply(res["//Hsp_align-len"], xmlValue),
		Hsp_qseq = sapply(res["//Hsp_qseq"], xmlValue),
		Hsp_hseq = sapply(res["//Hsp_hseq"], xmlValue)
		)))
	
	# in blast output, field named Hit_def contains a long string formatted <gene name><[organism name]>
	# we are interested in the gene name, so the below extracts from a string, everything before the first occurence of "["  
	annotation = str_split(as.character(blast_out$Hit_def[1]), fixed("["))[[1]][1]

	# save the blastp table as output file
	write.table(file=paste(outdir,paste(PC, gsub(' ','_',annotation),".txt", sep = "_"), sep = "/"),x= blast_out,quote=F,row.names=T,col.names=T,sep="\t",append=F)
	print('blast-out table done')

	# Here we create a fasta file from the blast output, extracting organism name from Hit_def (text between "[" and "]" )
	fasta_frame <- data.frame(paste(blast_out$Hit_accession,sub("\\].*", "", sub(".*\\[", "", blast_out$Hit_def)),sep='_'), blast_out$Hsp_hseq)
	df.fasta = dataframe2fas(fasta_frame,file=paste(outdir,paste(PC, annotation,".fasta", sep = "_"), sep = "/"))
	df.fasta = dataframe2fas(fasta_frame,file='df.fasta')
	print('to fasta done')

	# extract genus name from Hit_def, reformat it, make a dataframe from it and merge with TaxTable (ProkComp), ask for file if I forget! 
	closestHits = sub("\\].*", "", sub(".*\\[", "", blast_out$Hit_def))
	genera =  sub("\\ .*", "", closestHits)
	TaxonomyBlast = data.frame(cbind('hit'= closestHits,'Genus_name'= genera))
	TaxonomyBlast <- join(TaxonomyBlast, TaxTable, by = "Genus_name")
	TaxonomyBlast$hit2 = paste(blast_out$Hit_accession,TaxonomyBlast$hit,sep='_')
	subGeneInfo = data.frame(table(TaxonomyBlast $Genus_name),'gene'=paste(PC,annotation,sep='_'))
	
	#Store the taxonomical information to compare all your input sequences!!!
	OverallTaxonData= rbind(OverallTaxonData, subGeneInfo)
	topHit = c(topHit,as.character(TaxonomyBlast$hit[2]))
	
	#Closure of tryCatch({ which prevents ERROR's to stop the loop
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


topHit = as.data.frame(table(topHit))
rownames(topHit) = topHit$topHit
colnames(topHit) = c('topHit','query')

th = data.frame(topHit[,'query'])
rownames(th)= topHit$topHit
colnames(th) = c('query')

return(th)
}

#-------------------------------------------------------------------------------
# 	FROM ANVIO
#-------------------------------------------------------------------------------


blastpANVIO <- function(genes_of_interest, GenusOfInterest, outdir ){
OverallTaxonData = rbind()
topHit = c()

for(PC in genes_of_interest){
#for(seqnr in 1:length(InputFasta)){
	tryCatch({
	print(paste("Currently starting the analysis for:",PC,sep=' '))
	
	res = blastSequences(as.character(ANVIO[ANVIO$protein_cluster_id== PC,'aa_sequence'][1]), "nr", program="blastp", 	as="XML",hitListSize="3",timeout=2000)
		print('res done')
		
	# blastSequences returns output in XML format, lukily it is convenently structured 
	blast_out = data.frame(cbind(cbind(	
		Hit_num = sapply(res["//Hit_num"], xmlValue),
		Hit_id= sapply(res["//Hit_id"], xmlValue), 
		Hit_def = sapply(res["//Hit_def"], xmlValue),
		Hit_accession = sapply(res["//Hit_accession"], xmlValue),
		Hit_len= sapply(res["//Hit_len"], xmlValue)),
		cbind(	
		Hsp_bit_score= sapply(res["//Hsp_bit-score"], xmlValue), 
		Hsp_score = sapply(res["//Hsp_score"], xmlValue),
		Hsp_evalue = sapply(res["//Hsp_evalue"], xmlValue), 
		Hsp_identity = sapply(res["//Hsp_identity"], xmlValue),
		Hsp_gaps = sapply(res["//Hsp_gaps"], xmlValue), 
		Hit_len= sapply(res["//Hit_len"], xmlValue),
		Hsp_align_len = sapply(res["//Hsp_align-len"], xmlValue),
		Hsp_qseq = sapply(res["//Hsp_qseq"], xmlValue),
		Hsp_hseq = sapply(res["//Hsp_hseq"], xmlValue)
		)))
	
	# in blast output, field named Hit_def contains a long string formatted <gene name><[organism name]>
	# we are interested in the gene name, so the below extracts from a string, everything before the first occurence of "["  
	annotation = str_split(as.character(blast_out$Hit_def[1]), fixed("["))[[1]][1]

	# save the blastp table as output file
	write.table(file=paste(outdir,paste(PC, gsub(' ','_',annotation),".txt", sep = "_"), sep = "/"),x= blast_out,quote=F,row.names=T,col.names=T,sep="\t",append=F)
	print('blast-out table done')

	# Here we create a fasta file from the blast output, extracting organism name from Hit_def (text between "[" and "]" )
	fasta_frame <- data.frame(paste(blast_out$Hit_accession,sub("\\].*", "", sub(".*\\[", "", blast_out$Hit_def)),sep='_'), blast_out$Hsp_hseq)
	df.fasta = dataframe2fas(fasta_frame,file=paste(outdir,paste(PC, annotation,".fasta", sep = "_"), sep = "/"))
	df.fasta = dataframe2fas(fasta_frame,file='df.fasta')
	print('to fasta done')

	# extract genus name from Hit_def, reformat it, make a dataframe from it and merge with TaxTable (ProkComp), ask for file if I forget! 
	closestHits = sub("\\].*", "", sub(".*\\[", "", blast_out$Hit_def))
	genera =  sub("\\ .*", "", closestHits)
	TaxonomyBlast = data.frame(cbind('hit'= closestHits,'Genus_name'= genera))
	TaxonomyBlast <- join(TaxonomyBlast, TaxTable, by = "Genus_name")
	TaxonomyBlast$hit2 = paste(blast_out$Hit_accession,TaxonomyBlast$hit,sep='_')
	subGeneInfo = data.frame(table(TaxonomyBlast $Genus_name),'gene'=paste(PC,annotation,sep='_'))
	
	#Store the taxonomical information to compare all your input sequences!!!
	OverallTaxonData= rbind(OverallTaxonData, subGeneInfo)
	topHit = c(topHit,as.character(TaxonomyBlast$hit[2]))
	
	#Closure of tryCatch({ which prevents ERROR's to stop the loop
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


topHit = as.data.frame(table(topHit))
rownames(topHit) = topHit$topHit
colnames(topHit) = c('topHit','query')

th = data.frame(topHit[,'query'])
rownames(th)= topHit$topHit
colnames(th) = c('query')

return(th)
}






#-------------------------------------------------------------------------------
# 	FROM ANVIO
#-------------------------------------------------------------------------------

genes_of_interest = c('PC_00005278','PC_00005173','PC_00005166','PC_00005288','PC_00005191','PC_00005167','PC_00005226','PC_00005295','PC_00005307','PC_00005287','PC_00005346','PC_00005242','PC_00005368','PC_00005190','PC_00005168','PC_00005201','PC_00005351','PC_00005397','PC_00005284','PC_00004976','PC_00005029','PC_00005080','PC_00005088','PC_00004947')


OverallTaxonData2 = blastpANVIO(genes_of_interest, GenusOfInterest, outdir)





#-------------------------------------------------------------------------------
# 	FROM FASTA
#-------------------------------------------------------------------------------
#####

OverallTaxonData2 = blastpFASTA(InputFasta, GenusOfInterest, outdir)
chordDiagram(t(OverallTaxonData2))

topHit = as.data.frame(table(OverallTaxonData[,c(1)]))
rownames(topHit) = topHit$topHit
colnames(topHit) = c('topHit','query')

th = data.frame(topHit[,'query'])
rownames(th)= topHit$topHit
colnames(th) = c('query')

chordDiagram(t(th))

######
topHit = as.data.frame(table(OverallTaxonData[,c(1)]))
rownames(topHit) = topHit$topHit
colnames(topHit) = c('Genus_name','query')

th = data.frame(topHit)
th <- join(topHit, TaxTable, by = "Genus_name")
th = th[,c('Class','query')]

topHit = as.data.frame(table(th[,c(1)]))
rownames(topHit) = th$Var1
colnames(topHit) = c('topHit','query')

th = data.frame(topHit[,'query'])
rownames(th)= topHit$topHit
colnames(th) = c('query')


chordDiagram(t(th))

######

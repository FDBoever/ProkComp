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

SilvaTree = read.tree('LTPs128_SSU_tree.newick')


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


OverallTaxonData2 = blastpFASTA(InputFasta, GenusOfInterest, outdir)



blastpFASTA <- function(InputFasta, GenusOfInterest, outdir ){




OverallTaxonData = rbind()

#for(PC in genes_of_interest){
for(seqnr in 1:length(InputFasta)){
	PC = InputFasta@ranges@NAMES[seqnr]
	tryCatch({
	print(paste("Currently starting the analysis for:",PC,sep=' '))
	#res = blastSequences(as.character(ANVIO[ANVIO$protein_cluster_id== PC,'aa_sequence'][1]), "nr", program="blastp", 	as="XML",hitListSize="100",timeout=2000)
	
	
	#performing the blast search, generally slow, (tweak maxium waiting time by changing 'timeout' value )
	res = blastSequences(as.character(InputFasta[seqnr]), "nr", program="blastp", as="XML",hitListSize="100",timeout=2000)
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


	# reload the fastafile as Biostrings AAStringSet object
	sequences <- Biostrings::readAAStringSet('df.fasta',format='fasta')
	print('load sequences done')


	#align the sequences using CLUSTAL 2.1
	aln <- msa(sequences)
	export.fasta(msaConvert(aln, type=c("bios2mds::align")), outfile = "alignment.fa", ncol = 60, open = "w")
	print(aln, show="complete")
	aln <- msaConvert(aln, type="seqinr::alignment")
	print('align sequences done')


	# dirty but fast: distance matrix and NJ tree reconstruction, and assesment of monophyly of a genus (CHANGE MARINOBACTER); 
	dist.aln <- dist.alignment(aln, "identity")
	njTree <- njs(dist.aln )
	njTreeMR = ladderize(midpoint.root(njTree), right = FALSE)
	monophyl = is.monophyletic(phy = njTree, tips = njTreeMR$tip.label[grepl(GenusOfInterest,njTreeMR$tip.label)])
	print('tree building done')


	# OPTIONAL: Prints a tree and save as pdf; 
	pdf(paste(outdir,paste(PC, gsub(' ','_',annotation),monophyl,".pdf", sep = "_"), sep = "/"))
	plot(njTreeMR, font = 1, cex = 0.3,main=annotation)
	dev.off()
	print('file saving done')
	
	# extract genus name from Hit_def, reformat it, make a dataframe from it and merge with TaxTable (ProkComp), ask for file if I forget! 
	closestHits = sub("\\].*", "", sub(".*\\[", "", blast_out$Hit_def))
	genera =  sub("\\ .*", "", closestHits)
	TaxonomyBlast = data.frame(cbind('hit'= closestHits,'Genus_name'= genera))
	TaxonomyBlast <- join(TaxonomyBlast, TaxTable, by = "Genus_name")
	TaxonomyBlast$hit2 = paste(blast_out$Hit_accession,TaxonomyBlast$hit,sep='_')
	subGeneInfo = data.frame(table(TaxonomyBlast $Genus_name),'gene'=paste(PC,annotation,sep='_'))
	
	#format a list of genera containing the strain names, needed for groupOTU function in ggtree 
	l <- list();	
	for(genus in TaxonomyBlast$Genus_name){
		l[[genus]] = as.character(TaxonomyBlast[TaxonomyBlast$Genus_name==genus,'hit2'])
	}
	annotatedTree <- groupOTU(njTreeMR, l)
	#ggtree(annotatedTree, aes(color=group), layout='circular') + geom_tiplab(size=1, aes(angle=angle))

	genusTip = annotatedTree$tip.label[grepl(GenusOfInterest,annotatedTree$tip.label)]
	rca = getMRCA(annotatedTree,tip= genusTip)
	
	#Several ways of visualising the alignment and NJ tree
	p1 = msaplot(ggtree(njTreeMR)+geom_tiplab(size=2)+geom_treescale(offset=-1)+ggtitle(annotation), 'alignment.fa',offset=0.2,width = 2)
	p1.2 = msaplot(ggtree(njTreeMR)+geom_treescale(offset=-1)+ggtitle(annotation), 'alignment.fa',width = 4)
	p2 = ggplot() + geom_bar(aes(y = Freq, x = annotation, fill = Var1), data = data.frame(table(TaxonomyBlast $Genus_name)),stat="identity")
	p3 = ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,) +geom_treescale()+geom_hilight(node= rca, fill="gray")+geom_cladelabel(node= rca, label= GenusOfInterest,offset.text=0.05)+theme_tree()
	p4 = ggtree(annotatedTree)+geom_tiplab(aes(color= group),size=2)+geom_tippoint(aes(color= group), alpha=1) +geom_treescale() 
	p6= ggtree(annotatedTree)+geom_tippoint(aes(color= group), alpha=1) +geom_treescale() + xlim(0, 1) +geom_hilight(node= rca, fill="lightgray")+geom_cladelabel(node= rca, label= GenusOfInterest)+theme_tree()
	p3.1 = ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,)


	
	p3  = ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,) +geom_treescale()+geom_tiplab(aes(color= group,angle=angle),size=2,align=TRUE,linesize=.2)+geom_hilight(node= rca, fill="gray")

	#visualise and save the prefered graphic
	pdf(paste(outdir,paste('multi_',PC, gsub(' ','_',annotation),monophyl,".pdf", sep = "_"), sep = "/"))
	multiplot(p1.2,p4,p6,p3,cols=4)
	dev.off()
	
	#Store the taxonomical information to compare all your input sequences!!!
	OverallTaxonData= rbind(OverallTaxonData, subGeneInfo)
	
	d1 <- data.frame(id= annotatedTree $tip.label, val = as.numeric(as.character(blast_out$Hsp_bit_score)))
	d2 <- data.frame(id= annotatedTree $tip.label, value = as.numeric(as.character(blast_out$Hsp_evalue)))
	p = ggtree(annotatedTree)+geom_tippoint(aes(color= group), alpha=1) +geom_treescale()+geom_tiplab(aes(color= group),size=2)
	p2 <- facet_plot(p, panel="bit-score", data=d1, geom=geom_point, aes(x=val), color='firebrick')
	facet_plot(p2, panel='identity', data=d2, geom= geom_point, aes(x=value), color='steelblue') + theme_tree2()
	
	#blast_pca = blast_out[,c('Hit_len','Hsp_bit_score','Hsp_score','Hsp_evalue', 'Hsp_identity','Hsp_gaps','Hit_len.1','Hsp_align_len')]
	#indx <- sapply(blast_pca, is.factor)
	#blast_pca[indx] <- lapply(blast_pca[indx], function(x) as.numeric(as.character(x)))


	#res.pca <- PCA(blast_pca, graph = FALSE)
	#eig.val <- get_eigenvalue(res.pca)
	#fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
	#var <- get_pca_var(res.pca)

	#TaxOut = cbind(blast_pca,'Genus_name'=genera)
	#TaxOut = join(TaxOut, TaxTable, by = "Genus_name")

	#pdf(paste(outdir,paste('PCA_',gsub(' ','_',annotation),monophyl,".pdf", sep = "_"), sep = "/"))
	#fviz_pca_biplot(res.pca, 
                # Individuals
    #            geom.ind = "point",
    #            fill.ind = genera, col.ind = "white",
    #            pointshape = 21, pointsize = 2,
    #            addEllipses = TRUE,
    #            # Variables
    #            alpha.var ="contrib", col.var = "contrib",
    #            gradient.cols = "RdBu"
    #            )+ labs(fill = "Species", color = "Contrib", alpha = "Contrib") # Change legend title	
  	#dev.off()
	
	p3.1 = ggtree(annotatedTree,layout='circular')+geom_tippoint(aes(color= group), alpha=1,)
	colint = 200
	for(genus in unique(genera)){
	genusTip = annotatedTree$tip.label[grepl(paste(genus,' ',sep=''),annotatedTree$tip.label)]
	rca = getMRCA(annotatedTree,tip= genusTip)
	print(rca)
	p3.1 = p3.1 + geom_hilight(node= rca, fill=colors()[colint])+geom_cladelabel(node= rca, label= genus,offset.text=0.05)
	colint = colint + 10 
	}

	p6= ggtree(annotatedTree)+geom_tippoint(aes(color= group), alpha=1) +geom_treescale() + xlim(0, 1)+ theme_tree() + theme(legend.position="none")
	colint = 200
	for(genus in  unique(TaxTable$Genus_name[TaxTable$Genus_name %in% unique(genera)])[3:19]){
	genusTip = annotatedTree $tip.label[grepl(paste(genus,' ',sep=''), annotatedTree $tip.label)]
	rca = getMRCA(annotatedTree,tip= genusTip)
	print(rca)
	p6 = p6 + geom_hilight(node= rca, fill=colors()[colint])+geom_cladelabel(node= rca, label= genus,offset.text=0.05)
	colint = colint + 10 
	}

	
	pdf(paste(outdir,paste('circular_annotated_',gsub(' ','_',annotation),monophyl,".pdf", sep = "_"), sep = "/"))
	multiplot(p3.1,p6,cols=2)
	dev.off()
	
	silva =  ggtree(SilvaTree,layout='circular',size=0.3)
	colint = 200
	for(genus in unique(TaxTable$Genus_name[TaxTable$Genus_name %in% unique(genera)])){
	genusTip = SilvaTree $tip.label[grepl(paste(genus,'_',sep=''), SilvaTree $tip.label)]
	rca = getMRCA(SilvaTree,tip= genusTip)
	print(rca)
	silva = silva + geom_hilight(node= rca, fill=colors()[colint])+geom_cladelabel(node= rca, label= genus,offset.text=0.05)
	colint = colint + 10 
	}
	
	pdf(paste(outdir,paste('SILVA',gsub(' ','_',annotation),monophyl,".pdf", sep = "_"), sep = "/"))
	multiplot(silva)
	dev.off()
	
	#Closure of tryCatch({ which prevents ERROR's to stop the loop
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
return(OverallTaxonData)
}


# compare different taxonomic levels 
#ggplot() + geom_bar(position = "fill",stat = "identity",aes(y = Freq, x = sort(Var1, Freq)), data = OverallTaxonData,stat="identity")+ guides(fill=FALSE)+coord_flip()
#ggplot() + geom_bar(aes(y = Freq, x = gene, fill = Var1), data = OverallTaxonData,stat="identity")+ guides(fill=FALSE)+coord_flip()

colnames(OverallTaxonData2) =  c('Genus_name','Freq','gene')

OverallTaxonData = join(OverallTaxonData2, TaxTable, by = "Genus_name")
q1 = ggplot() + geom_bar(position = "fill",aes(y = Freq, x = gene, fill = Order), data = OverallTaxonData,stat="identity")+ theme(axis.text.y = element_text(size=5))+coord_flip()+scale_fill_viridis(discrete=TRUE)
q2 = ggplot() + geom_bar(position = "fill",aes(y = Freq, x = gene, fill = Class), data = OverallTaxonData,stat="identity")+theme(axis.text.y = element_text(size=5))+coord_flip()+scale_fill_viridis(discrete=TRUE)
q3 = ggplot() + geom_bar(position = "fill",aes(y = Freq, x = gene, fill = Phylum), data = OverallTaxonData,stat="identity")+theme(axis.text.y = element_text(size=5))+coord_flip()+scale_fill_viridis(discrete=TRUE)


multiplot(q1,q2,q3)


ggplot() + geom_bar(position = "fill",aes(y = Freq, x = gene, fill = Order), data = OverallTaxonData,stat="identity")+coord_flip()


#-------------- blast_PCA -----------------

library("FactoMineR")
library("factoextra")
library(reshape)
library(heatmap.2)
library("corrplot")


silva =  ggtree(SilvaTree,layout='circular',size=0.3)

colint = 200
for(genus in unique(TaxTable$Genus_name[TaxTable$Genus_name %in% unique(genera)])){
	genusTip = SilvaTree $tip.label[grepl(paste(genus,'_',sep=''), SilvaTree $tip.label)]
	rca = getMRCA(SilvaTree,tip= genusTip)
	print(rca)
	silva = silva + geom_hilight(node= rca, fill=colors()[colint])+geom_cladelabel(node= rca, label= genus,offset.text=0.05)
	colint = colint + 10 
}






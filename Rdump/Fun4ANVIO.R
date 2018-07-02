




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


####################################################################################


FindSeqInAnvio = function(sequence, ANVIOdb){
	#remove the dashed introduced by algining
	SelectedANVIO = ANVIOdb[grepl(sequence, ANVIOdb aa_seq_trim,ignore.case=TRUE),]
	
	print(cat(':: Query Sequence ::\n',sequence,'\n- Genome :', as.character(SelectedANVIO$genome_name),'\n- PC :', as.character( SelectedANVIO$protein_cluster_id),'\n- COG_CAT :', as.character( SelectedANVIO$COG_CATEGORY_ACC),'\n- COG_FUNC :', as.character( SelectedANVIO$COG_FUNCTION_ACC),'\n - Annotation :', as.character( SelectedANVIO$COG_FUNCTION)))	
	return(as.character(SelectedANVIO$protein_cluster_id))
}

FindPCinPAN = function(PC,panmatrix){
	return(panmatrix[grepl(paste(PC,collapse='|'),rownames(panmatrix)),])
}

FindPCinANVIO= function(PC,ANVIO){
	return(ANVIO[grepl(paste(PC,collapse='|'),ANVIO$protein_cluster_id),])
}

FindFuncinANVIO= function(Func,ANVIO){
	return(ANVIO[grepl(paste(Func,collapse='|'), ANVIO$COG_FUNCTION,ignore.case=TRUE),])
}

Fasta2Pan = function(InputFasta,AnvioData,ANVIO){
	PC.tracker = c()
	for(seqnr in 1:length(InputFasta)){
		name = InputFasta@ranges@NAMES[seqnr]
		print(paste("Currently starting the analysis for:", name,sep=' '))	
		testPC = FindSeqInAnvio(as.character(InputFasta[seqnr]), ANVIO)
		PC.tracker = c(PC.tracker,testPC)
	}
	testPan = FindPCinPAN(PC.tracker, AnvioData)
	return(testPan)
}

#UniprotHeader = funciton(InputFasta){
#	InputFasta@ranges@NAMES
#}

PC2Tree= function(PCs,ANVIO,outdir){
	for(PC in PCs){
		SelectedANNIO = 	ANVIO[grepl(PC,ANVIO$protein_cluster_id),]
		fasta_frame <- data.frame(SelectedANNIO $genome_name, SelectedANNIO $aa_sequence)

		df.fasta = dataframe2fas(fasta_frame,file=paste(outdir,paste(PC,".fasta", sep = "_"), sep = "/"))
		df.fasta = dataframe2fas(fasta_frame,file='df.fasta')
		print('to fasta done')

		sequences <- Biostrings::readAAStringSet('df.fasta',format='fasta')
		print('load sequences done')

		#align the sequences using CLUSTAL 2.1
		aln <- msa(sequences)
		aln <- msaConvert(aln, type="seqinr::alignment")
		print('align sequences done')


		# dirty but fast: distance matrix and NJ tree reconstruction, and assesment of monophyly of a genus (CHANGE MARINOBACTER); 
		dist.aln <- dist.alignment(aln, "identity")
		njTree <- njs(dist.aln )
		njTree$edge.length = njTree$edge.length + 0.0000001
		njTree = drop.tip(njTree,njTree$tip.label[!grepl(paste(rownames(CheckM3_annotated),collapse='|'),njTree$tip.label)])

		p <- ggtree(midpoint.root(njTree),ladderize=TRUE)+geom_treescale() + xlim(0,0.5)
			p1 <- p %<+% CheckM3_annotated[,c(1:10)] + geom_tippoint(aes(color=group2))+ geom_tiplab(aes(color=group2))+scale_colour_manual(values=colors)
		print(p1)
	}
	#return(ANVIO[grepl(paste(PC,collapse='|'),ANVIO$protein_cluster_id),])
}

GrepAnvioHEAT = function(toSearch){
	testPan = FindPCinPAN(as.character(unique(FindFuncinANVIO(paste(toSearch,collapse='|'), ANVIO)[,'protein_cluster_id'])), AnvioData)
	testANVIO = FindPCinANVIO(rownames(testPan), ANVIO)
	dedtestANVIO = testANVIO[,c('protein_cluster_id','COG_CATEGORY_ACC','COG_FUNCTION_ACC', 'COG_FUNCTION')]
	dedtestANVIO  = dedtestANVIO[!duplicated(dedtestANVIO$protein_cluster_id), ]
	rownames(dedtestANVIO)= dedtestANVIO$protein_cluster_id
	dedtestANVIO = dedtestANVIO[rownames(testPan),]
	rownames(testPan) = paste(dedtestANVIO$protein_cluster_id, dedtestANVIO$COG_FUNCTION_ACC, dedtestANVIO$COG_FUNCTION,sep=' - ')
	print(heatmap.2(as.matrix((testPan)),trace='none',margins=c(20,20),cexRow=0.5,cexCol =0.5,col=colorRampPalette(c("white", "darkblue", "red",'green'))(n = 20)))	
}



FASTA2Tree= function(path){
		sequences <- Biostrings::readAAStringSet(path,format='fasta')
		print('load sequences done')

		#align the sequences using CLUSTAL 2.1
		aln <- msa(sequences)
		aln <- msaConvert(aln, type="seqinr::alignment")
		print('align sequences done')


		# dirty but fast: distance matrix and NJ tree reconstruction, and assesment of monophyly of a genus (CHANGE MARINOBACTER); 
		dist.aln <- dist.alignment(aln, "identity")
		njTree <- njs(dist.aln )
		njTree$edge.length = njTree$edge.length + 0.0000001

		p <- ggtree(midpoint.root(njTree),ladderize=TRUE)+geom_treescale() + xlim(0,2)
			p1 <- p + geom_tippoint()+ geom_tiplab(size=2)
		print(p1)

	#return(ANVIO[grepl(paste(PC,collapse='|'),ANVIO$protein_cluster_id),])
}



####################################################################################


#Visialise a NJ tree on the basis of an input fasta file
FASTA2Tree('~/Genomics/alkB_all.fasta')
FASTA2Tree('~/Genomics/alkB_all_2.fasta')




ANVIO$aa_seq_trim = gsub('-','', ANVIO$aa_seq)

InputFasta = readAAStringSet('~/Genomics/Glycolate_fasta.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/Glycolate_HP15_2.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/full_glycolate.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/TRAP_mannitol.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/alkB_all.fasta',format='fasta')

~/Genomics/catechol_ortho.fasta
~/Genomics/catechol_meta.fasta
~/Genomics/full_denitrification.fasta


InputFasta = readAAStringSet('~/Genomics/glycogen_hp15.fasta',format='fasta')


SeqOfInterest = 'MADLSQKLKEQVLQARKTGQKLNIEGGGTKAFMGRAADTDAGSLKMGEHTGIVEYHPVELILTARAGTTLAEIEATLAEEGQCLHFEPPRFGDGSTLGGTLACNLSGPARPWSGSVRDQVLGIRLLNGKGEHLRFGGQVMKNVAGYDTSRLQAGAMGTLGAITEISLKVMPKPAMSLTLVQDMAMDKVIHYMNSRAKEPKPITGAAWVDGKVYLRLAGAKSAVEATAEKWSGEVMEQGDEFWQQLQHMQHDFFAGDEPLWRFSVGSTAETPDLDGPWLIDWAGSQRWYRGEADIAQLETMAQRAGGQVSLFRGGDRTDEVMHHQPRALKEIQQRLKKSFDPDGIFNPGRLYSWL'

SeqOfInterest =  'MVQRISQAFLVMFVISVIAFAIQDGLGDPLQEMVGMSVSEEEREAIREELGLNDPMVVQYLRFAGNALQGDLGNSYFYSKPTLEVIAEHLPATLELVIGASLIIVFLSVPIGVYAAIRPQAWLSKFFMGVSTVGISIPVFLTAIVLIQLFSIGVTVNWFPSDTGWGQWLNDFFSTEGGLPSYGRGDELTHLFGTWESGFFTGTGLMHLVLPSVSLASIMLPLFIRLIRAEMMEVLQSDYVRYARAKGLSAGRINFLHALKNTMLPVITVGGVQIGIMVAYTILTETVFQWPGVGLMFLEAITRSDIPLIVAYLMVVGLIFVVTNTLVDLIYGLVNPTVKLTGKKA'



testPan = FindPCinPAN(testPC, AnvioData)
testPan = FindPCinPAN(names(tail(sort(rowSums(AnvioData)),5)), AnvioData)


#Grep the pangenome file for function(s) and visualise a heatmap
GrepAnvioHEAT('lactate')
GrepAnvioHEAT('flag')
GrepAnvioHEAT('lactate')
GrepAnvioHEAT('biotin')
GrepAnvioHEAT('efflux pump')
GrepAnvioHEAT('desaturase')
GrepAnvioHEAT('mannitol')
GrepAnvioHEAT('mannose')
GrepAnvioHEAT('oxidored')
GrepAnvioHEAT(c('precorrin','cobalamin','cobyrnic'))
GrepAnvioHEAT('sulfate')
GrepAnvioHEAT('copper')
GrepAnvioHEAT('efflux system')
GrepAnvioHEAT('Fatty Acid desaturase')
GrepAnvioHEAT('glyoxylate')




InputFasta = readAAStringSet('~/Genomics/urea_MAQU.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/HP15_flagellar_operon.org',format='fasta')
InputFasta = readAAStringSet('~/Genomics/oligopeptide_hp15.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/catechol_ortho.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/catechol_meta.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/full_denitrification.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/all_alkB.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/xantho.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/glcD_all.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/glcE_all.fasta',format='fasta')
InputFasta = readAAStringSet('~/Genomics/glcF_all.org',format='fasta')
InputFasta = readAAStringSet('~/Genomics/nuo_blaat.fasta',format='fasta')


testFasta2Pan = Fasta2Pan(InputFasta, AnvioData, ANVIO)
testANVIO = FindPCinANVIO(rownames(testFasta2Pan), ANVIO)

dedtestANVIO = testANVIO[,c('protein_cluster_id','COG_CATEGORY_ACC','COG_FUNCTION_ACC', 'COG_FUNCTION')]
dedtestANVIO  = dedtestANVIO[!duplicated(dedtestANVIO$protein_cluster_id), ]
rownames(dedtestANVIO)= dedtestANVIO$protein_cluster_id
dedtestANVIO = dedtestANVIO[rownames(testFasta2Pan),]

rownames(testFasta2Pan) = paste(dedtestANVIO$protein_cluster_id, dedtestANVIO$COG_FUNCTION_ACC, dedtestANVIO$COG_FUNCTION,sep=' - ')

heatmap.2(as.matrix((testFasta2Pan)),trace='none',margins=c(30,30),cexRow=0.5,cexCol =0.5)	



as.character(unique(FindFuncinANVIO('lactate', ANVIO)[,'protein_cluster_id']))


PC2Tree(c('PC_00000494','PC_00003588'),ANVIO, '~/')




		sequences <- Biostrings::readAAStringSet(path,format='fasta')
		print('load sequences done')

		#align the sequences using CLUSTAL 2.1
		aln <- msa(sequences)
		aln <- msaConvert(aln, type="seqinr::alignment")
		print('align sequences done')


# dirty but fast: distance matrix and NJ tree reconstruction, and assesment of monophyly of a genus (CHANGE MARINOBACTER); 
dist.aln <- dist.alignment(aln, "identity")
njTree <- njs(dist.aln )
njTree$edge.length = njTree$edge.length + 0.0000001

p <- ggtree(midpoint.root(njTree),ladderize=TRUE)+geom_treescale() + xlim(0,2)
p1 <- p + geom_tippoint()+ geom_tiplab(size=2)
	print(p1)

mnjTree = midpoint.root(njTree)
mnjTree$tip.label[grep('Marinobacter ',mnjTree$tip.label)]


njTree

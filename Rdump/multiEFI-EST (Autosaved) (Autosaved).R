#	EFI-EST
#####################

# Rationale
#	based on rudimental phylogeny and explorations on uniprot and such, we identified the high diversity of alkB genes ~ several clades of alkB "subfamilies"
#	we descided to use EFI-EST to visualise potential subfamilies and closely related sequences for each of the alkb's in marinobacter
#		for this we:
#			1) we ran EFI-EST "sequence BLAST" on each of the following aa sequences consecutatively
#					MAQU_0440, MAQU_0610, MAQU_2843, MDG893_03170, ABA45_07610, ABO_0122, ABO_2707
#			2) in cytoscape, for each of the networks, we exported the "node table" to csv files
#			3) imported the csv files into R and export all the unique accession numbers to a local txt file
#			4) this file was then uploaded onto EFI-EST, and the final network was calculated online

####################

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

####################



# Load the csv files (EFI-EST generated)
#AlkB genes
maqu_0440 = read.csv('~/Downloads/0440.csv')
maqu_2843 = read.csv('~/Downloads/MAQU_2843_EFI-EST.csv')
maqu_0610 = read.csv('~/Downloads/maqu_0610_EFI-EST.csv')
MDG893_03170 = read.csv('~/Downloads/DG.csv')
ABA45_07610 = read.csv('~/Downloads/ABA45_07610_EFI-est.csv')
ABO_2707 = read.csv('~/Downloads/ABO_2707_EFI-EST.csv')
ABO_0122 = read.csv('~/Downloads/ABO_0122_EFI_EST.csv')

#ATPsynthase in NKSG1
MSNKSG1_00556 = read.csv('~/Downloads/MSNKSG1_00556.csv')
MSNKSG1_14497 = read.csv('~/Downloads/MSNKSG1_14497.csv')
MSNKSG1_04306 = read.csv('~/Downloads/MSNKSG1_04306.csv')


# ALKB - select the unique ACCESSION numbers from all the files combined
selectedACCession = unique(c(as.character(maqu_0440 $name),as.character(maqu_2843 $name),as.character(maqu_0610 $name),as.character(MDG893_03170 $name),as.character(ABA45_07610 $name),as.character(ABO_2707 $name),as.character(ABO_0122 $name)))

# ATPASE - select the unique ACCESSION numbers from all the files combined
selectedACCessionAPTase = unique(c(as.character(MSNKSG1_00556 $name),as.character(MSNKSG1_14497 $name),as.character(MSNKSG1_04306 $name)))





#save to file (which is ready to be uplaoded to EFI-EST)
write(selectedACCession,'data.txt')
write(selectedACCessionAPTase,'dataATP.txt')

selectCols = c("SUID","BRENDA.ID","CAZY.Name","Class","Description","EC","Family","Gene.Name","Genus","GO.Term","HMP.Body.Site","HMP.Oxygen","IPRO","KEGG.ID","Kingdom","name","NCBI.IDs","Order","Organism","P01.gDNA","selected","Sequence.Length","shared.name","Species","STRING.ID","Superkingdom","Swissprot.Description","Taxonomy.ID","UniProt.Annotation.Status")


selectedTaxonIDs = uberEST$Taxonomy.ID
write(selectedTaxonIDs,'dataTaxonID.txt')

~/Downloads/dbfetch_res1.fasta

# TO DOWNLOAD THE RALTED AA SEQUENCES, (rather creative way of doing this), wget the below
#paste('https://www.uniprot.org/uniprot/?query=',paste(selectedACCession,collapse='+or+'),'&format=fasta ',sep='')
#paste('https://www.uniprot.org/uniprot/?query=',paste(selectedACCession[1:10],collapse='+or+'),'&format=fasta ',sep='')


####################
####################
####################


atpD_cytoscape = read.csv('~/Downloads/atpD_180_combined.csv')
alkB_cytoscape = read.csv('~/Downloads/alkB_full_combined.csv')

selectedACCession = atpD_cytoscape$name
selectedACCession = alkB_cytoscape$name

#selectedACCessionUNIPROT = paste('UNIPROT:',selectedACCession,sep='')
seqPerFasta = 200
nACC = length(selectedACCession)
nCycles = ceiling(nACC/seqPerFasta)

#UNIPROT is restricted in using 200 query items at the time, so we need some steps to get arount this

for(i in c(1:nCycles)){
	if(i!=nCycles){
		start=i* seqPerFasta-seqPerFasta +1
		end=i* seqPerFasta +1
	}else{
		start= i* seqPerFasta-seqPerFasta +1
		end= nACC
	}
	#print(selectedACCession[start:end])
	#download.file(paste('https://www.uniprot.org/uniprot/?query=',paste(selectedACCession[start:end],collapse='+or+'),'&format=fasta',sep=''),destfile=paste('~/part',i,'.fasta',sep=''))
	
	#print(selectedACCessionUNIPROT[start:end])

	print(paste('http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=',URLencode(paste(selectedACCession[start:end],collapse='%2C')),'&format=fasta&style=raw&Retrieve=Retrieve',sep=''))
	download.file(paste('http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=',URLencode(paste(selectedACCession[start:end],collapse='%2C')),'&format=fasta&style=raw&Retrieve=Retrieve',sep=''), paste('~/',i,'_uniprot_atpd.fasta',sep=''))

}



##### http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/

#### 
paste('http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=',paste(selectedACCession[start:end],collapse='%'),'&format=fasta&style=raw&Retrieve=Retrieve',sep='')




#write(c(paste('UNIPROT:', selectedACCession[start],'\n',sep=''),paste(selectedACCession[(start+1):end],collapse='\nUNPROT:')),'TESTTEST.txt')
#	print(c('UNIPROT:',paste(selectedACCession[start:end],collapse='\nUNPROT:')))




###########################

sequences <- Biostrings::readAAStringSet('~/1_uniprot_atpd.fasta',format='fasta')

sequence@ranges@NAMES = sub("\\_.*", "", sequence@ranges@NAMES)
sequence@ranges@NAMES = sub("UNIPROT:", "", sequence@ranges@NAMES)


path = "~/Documents/My Data/BRAZIL/Elections/"
out.file<-""
file.names <- dir(path, pattern =".txt")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i],header=TRUE, sep=";", stringsAsFactors=FALSE)
  out.file <- rbind(out.file, file)
}
 write.table(out.file, file = "cand_Brazil.txt",sep=";", 
             row.names = FALSE, qmethod = "double",fileEncoding="windows-1252")



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

###########################
alkB_cytoscape = read.csv('~/Downloads/alkB_full_combined.csv',stringsAsFactors=FALSE)


alkB_tree = read.tree('~/Genomics/alkB_fasttree.tre')
atpD_tree =read.tree('~/Genomics/atpD_fasttree.tre')  

alkB_tree=midpoint.root(alkB_tree)
atpD_tree =midpoint.root(atpD_tree)

alkB_tree$tip.label = sub("\\_.*", "", alkB_tree$tip.label)
alkB_tree$tip.label = sub("tr|", "", alkB_tree$tip.label)
alkB_tree$tip.label  = unique(unlist(strsplit(alkB_tree$tip.label, "[| ]")))

rownames(alkB_cytoscape) = alkB_cytoscape$name
ralkB = alkB_cytoscape
ralkB$name = as.character(ralkB$name)
#ralkB = ralkB[as.character(alkB_tree$tip.label), ]
#ralkB = ralkB[grepl(paste(alkB_tree$tip.label,collapse="|"),ralkB$name), ]


ralkB = alkB_cytoscape[grepl(paste(alkB_tree $tip.label,collapse="|"), alkB_cytoscape $name), ]
rownames(ralkB) = ralkB$name
ralkB = ralkB[alkB_tree $tip.label,]


groupInfo <- split(as.character(ralkB$name), ralkB$Genus)
chiroptera <- groupOTU(chiroptera, groupInfo)





p <- ggtree(alkB_tree,ladderize=TRUE,layout='circular')+geom_treescale() + xlim(0,2)

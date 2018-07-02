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




# Load the csv files (EFI-EST generated)
maqu_0440 = read.csv('~/Downloads/0440.csv')
maqu_2843 = read.csv('~/Downloads/MAQU_2843_EFI-EST.csv')
maqu_0610 = read.csv('~/Downloads/maqu_0610_EFI-EST.csv')
MDG893_03170 = read.csv('~/Downloads/DG.csv')
ABA45_07610 = read.csv('~/Downloads/ABA45_07610_EFI-est.csv')
ABO_2707 = read.csv('~/Downloads/ABO_2707_EFI-EST.csv')
ABO_0122 = read.csv('~/Downloads/ABO_0122_EFI_EST.csv')


# select the unique ACCESSION numbers from all the files combined
selectedACCession = unique(c(as.character(maqu_0440 $name),as.character(maqu_2843 $name),as.character(maqu_0610 $name),as.character(MDG893_03170 $name),as.character(ABA45_07610 $name),as.character(ABO_2707 $name),as.character(ABO_0122 $name)))

#save to file (which is ready to be uplaoded to EFI-EST)
write(selectedACCession,'data.txt')

selectCols = c("SUID","BRENDA.ID","CAZY.Name","Class","Description","EC","Family","Gene.Name","Genus","GO.Term","HMP.Body.Site","HMP.Oxygen","IPRO","KEGG.ID","Kingdom","name","NCBI.IDs","Order","Organism","P01.gDNA","selected","Sequence.Length","shared.name","Species","STRING.ID","Superkingdom","Swissprot.Description","Taxonomy.ID","UniProt.Annotation.Status")

uberEST = rbind(maqu_0440[, selectCols], mqu[, selectCols], alkB_DG[, selectCols], maqu_2843[, selectCols])
head(uberEST)



# TO DOWNLOAD THE RALTED AA SEQUENCES, (rather creative way of doing this), wget the below
paste('https://www.uniprot.org/uniprot/?query=',paste(selectedACCession,collapse='+or+'),'&format=fasta ',sep='')

paste('https://www.uniprot.org/uniprot/?query=',paste(selectedACCession[1:10],collapse='+or+'),'&format=fasta ',sep='')


nACC = length(selectedACCession)
nCycles = ceiling(nACC/200)

#UNIPROT is restricted in using 200 query items at the time, so we need some steps to get arount this

for(i in c(1:nCycles)){
	if(i!=nCycles){
		start=i*200-200+1
		end=i*200+1
	}else{
		start= i*200-200+1
		end= nACC
	}
	print(selectedACCession[start:end])
	download.file(paste('https://www.uniprot.org/uniprot/?query=',paste(selectedACCession[start:end],collapse='+or+'),'&format=fasta',sep=''),destfile=paste('~/part',i,'.fasta',sep=''))
}



download.file('https://www.uniprot.org/uniprot/?query=A0A1I4KFD0+or+A1TXS2+or+A0A137S556+or+A0A0K1UGC7+or+A0A0T6W2K0+or+A0A142FWT3+or+A0A2G1UGY0+or+A0A2N1X640+or+H8WCU7+or+A0A1H6ANQ2&format=fasta',destfile='~/part1.fasta')



dim(uberEST)
unique(uberEST$Description)
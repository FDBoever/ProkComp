############################################################################################################
#
# 16S diversity estimator
#
# -- RATIONALE
#	16S sequence alignments can be downloaded from the Silva-Webservice, I highly advice to mannulay check whether their are not duplicate strains and whether the alignment looks reasonable
#	Calculate pairwise sequence similarties using parSeqSim() from the protR package, as a measure for sequence diversity
# 	Load in other alignments, and calculte the same measure to have an idea how much of the 16S diversity is covered.
#
# --  KEY FUNCTIONS
# 	parSeqSim() - Parallelized pairwise similarity calculation with a list of protein sequences
#
#
# -- MARINOBACTER - A CASE STUDY 
#	The below example works with Marinobacter genus, in the framework of Comparative genomics
#	First we downloaded the Type strains from Silva, and look for 16S divergence measures
# 	we are interested to compare the measures against those obtained from our set of genomes (part of the comparison)
#
#
#
############################################################################################################



library('Biostrings')
library('gplots')
library('RColorBrewer')
library('ggfortify')

fasta = readDNAStringSet('~/Downloads/Silva_Marinobacter_T.fasta')

s1 <- fasta[1]
s2 <- fasta[2]

palign1 <- pairwiseAlignment(s1, s2)
palign1
pid(palign1)


#calculate all versus all sequence percent identities


row.tracker = c()
output.tracker = c()
print("... Starting calculations of  the percent sequence identity for each pairwise sequence alignment. ")
print("... May take a while depending on the number of sequences to compare")

for(i in 1:length(fasta)){
	s1 <- fasta[i]
	row.tracker = c()
	for(o in 1: length(fasta)){
		
	s2 <- fasta[o]
	palign1 <- pairwiseAlignment(s1, s2)
	pid.val = pid(palign1)
	row.tracker = c(row.tracker, pid.val)
	print(paste(i , "vs" , o , ":" , round(pid.val,4) , "% idenity"))
	}
	output.tracker = cbind(output.tracker, row.tracker)
}


colnames(output.tracker) = fasta@ranges@NAMES
rownames(output.tracker) = fasta@ranges@NAMES

heatmap.2(output.tracker,trace='none',col=colorRampPalette(brewer.pal(9, "RdBu"))(100))




#in case of the testing data, we observed a off-chart bad sequence, which indeed corresponds to a too small 16S fragment
# we kick him out as follows
output.tracker = output.tracker[!grepl('DSM 16072',rownames(output.tracker)),!grepl('DSM 16072',colnames(output.tracker))]







#in case of Silva derived fasta files, one could produce a dataframe looking like this:
DFnames = t(as.data.frame(strsplit(colnames(output.tracker),'Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Marinobacteraceae;Marinobacter;')))
colnames(DFnames) = c('accesion','name')

DFnames = data.frame(DFnames)
colnames(output.tracker) = paste(DFnames$name,DFnames$accesion)
rownames(output.tracker) = paste(DFnames$name,DFnames$accesion)



heatmap.2(output.tracker,trace='none',col=colorRampPalette(brewer.pal(9, "RdBu"))(100))

autoplot(prcomp(output.tracker))
autoplot(prcomp(output.tracker), label = TRUE, label.size = 3)


as.data.frame(strsplit(fasta@ranges@NAMES," "))
#	SelectNodes.R
# Script designed to throw out unwanted sequences from a fasta file
# Currently Optimised to use reference tree, and only keep the sequences which correspond to the tip.labels of the reference tree 

library(ape)
library(seqinr)
library(phytools)


referenceTree = read.tree('~/DATA/MarinobacterGenomics/2018_ProkComp/SCO_kde.fas.treefile')
RetainNodes = referenceTree$tip.label

alignmentsDIR = '~/DATA/MarinobacterGenomics/2018_ProkComp/alignments'
file_list <- list.files(alignmentsDIR,full.names=T)

for(alignment in file_list){
	sequences = read.fasta(alignment,seqtype='AA')
	trimmed_seq = sequences[grepl(paste(RetainNodes,collapse='|'), names(sequences))]
	write.fasta(trimmed_seq, names(trimmed_seq), file.out=paste(alignment,'trimmed.fa',sep=''))
}
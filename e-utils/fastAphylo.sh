#!/bin/sh

#	
#  	Java SE Development Kit 10.0.1
#


#	PASTE:		Combine into one; paste several fasta files underneath eachother
cat *.fasta >> aupA_combined.fasta

#	REMOVE:		Remove unwanted text, in this case "UNIPROT:"
sed -i -e 's/UNIPROT\://g' ./aupA_combined.fasta 

#	MUSCLE:	Run default MUSCLE ALIGNMENT using muscle 3.8.31
/Users/frederik/Genomics/muscle3.8.31_i86darwin64 -in aupA_combined.fasta -out aupB_aligned.fasta

#	BMGE:		TRIM alignment positions using BMGE 1.12
java -jar /Users/frederik/Genomics/BMGE-1.12/BMGE.jar -t AA -g 0.3 -i aupA_aligned.fasta -of aupA_aligned.fasta.bmge

#	FastTree:	Infer phylogeny using FastTree
/Users/frederik/Genomics/FastTree aupA_aligned.fasta.bmge > FastTree_aupA
/Users/frederik/Genomics/FastTree aupB_aligned.fasta.bmge > FastTree_aupB






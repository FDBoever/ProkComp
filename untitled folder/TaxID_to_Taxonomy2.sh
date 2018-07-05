#!/bin/bash

#-----------------------------------------------------------------------------------------
#	Taxon ID to taxonomy
#-----------------------------------------------------------------------------------------
#	Script will ute efetch to obtain taxonomy information of given NCBI TaxIDs
#
#	INPUT
#		<input file> 		a textfile with TaxIDs on single lines 
#	 
#	OUTPUT
#	Taxonomy_from_TaxID.txt a comma seperated text file with all TaxIDs and NCBI taxonomy
#
#	TO RUN:
# 	bash TaxID_to_Taxonomy.sh <input file>
#-----------------------------------------------------------------------------------------




#set -x 

TAXON_ID=$1
#mkdir ./out

for i in $(cat $TAXON_ID)
	do
		#echo $i
		TaxName = $(efetch -db taxonomy -id " " ${id} " " -format xml | xtract -pattern Taxon -first ScientificName)
		echo "$TaxName"
		#esearch -db nuccore -query "(("$TaxName "[Organism] AND 500:2000[Sequence Length])) AND 16S " | efetch -format fasta

		
		#outname=$(echo $i | sed 's/\//_/g' )
		#echo "### $i $outname ###" 
		#esearch -db nucleotide -query " CCAP AND "$i" " | efetch -format fasta >> ./out/CCAP_"$outname".fasta 
		#if [ -z "$outname" ]
		#	then
    #			rm ./out/CCAP_"$outname".fasta
	#	else
   # 			count=$(grep -c '>' ./out/CCAP_"$outname".fasta)
   # 			echo "CCAP $i has $count fasta record(s)"
#		fi
	done
exit 0


		#esearch -db nuccore -query "((Marinobacter sp. CP1[Organism] AND 500:2000[Sequence Length])) AND 16S " | efetch -format fasta



#efetch -db taxonomy -id 1761792,351348,1306954,1671721,1119533,1749259 -format xml | xtract -pattern Taxon -tab "," -first TaxId ScientificName \
#-group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" \
#-block "*/Taxon" -match "Rank:kingdom" -KING ScientificName \
#-block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName \
#-block "*/Taxon" -match "Rank:class" -CLSS ScientificName \
#-block "*/Taxon" -match "Rank:order" -ORDR ScientificName \
#-block "*/Taxon" -match "Rank:family" -FMLY ScientificName \
#-block "*/Taxon" -match "Rank:genus" -GNUS ScientificName \
#-group Taxon -tab "," -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&GNUS" >> Taxonomy_from_TaxID.txt
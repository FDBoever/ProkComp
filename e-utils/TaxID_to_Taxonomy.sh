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


TAXON_ID=$1
TAXON_IDs= cat $TAXON_ID | tr '\n' ','
echo "$TAXON_IDs"
echo "$TAXON_IDs"
echo "$TAXON_IDs"
echo "$TAXON_IDs"
echo $TAXON_ID | tr '\n' ','
echo $TAXON_ID | tr '\n' ','


efetch -db taxonomy -id echo "$TAXON_IDs" ' ' -format xml | xtract -pattern Taxon -tab "," -first TaxId ScientificName \
-group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" \
-block "*/Taxon" -match "Rank:kingdom" -KING ScientificName \
-block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName \
-block "*/Taxon" -match "Rank:class" -CLSS ScientificName \
-block "*/Taxon" -match "Rank:order" -ORDR ScientificName \
-block "*/Taxon" -match "Rank:family" -FMLY ScientificName \
-block "*/Taxon" -match "Rank:genus" -GNUS ScientificName \
-group Taxon -tab "," -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&GNUS" >> Taxonomy_from_TaxID.txt





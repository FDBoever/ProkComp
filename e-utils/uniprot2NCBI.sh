#!/bin/sh

echo '===================================================';
echo 'uniprot2NCBI.sh version beta 0.1';
echo '===================================================';

FILENAME=$1
echo $FILENAME

FILE=`basename ${FILENAME} .txt`
mkdir ${FILE};
echo '['`date "+%H:%M:%S"`'] - Loading File';

cat ${FILENAME} | while read line
	do
		echo '-------------------------------------------';
		echo "Processing $line"
		echo '-------------------------------------------';

		wget -qO- http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprot/$line/ > ${FILE}/${line}.out.txt
		NCBI=`cat ${FILE}/${line}.out.txt | grep RefSeq`
		EMBL=`cat ${FILE}/${line}.out.txt | grep EMBL`
		KEGG=`cat ${FILE}/${line}.out.txt | grep KEGG`
		PATRIC=`cat ${FILE}/${line}.out.txt | grep PATRIC`
		KO=`cat ${FILE}/${line}.out.txt | grep KO`
		PFAM=`cat ${FILE}/${line}.out.txt | grep Pfam`
		SubName=`cat ${FILE}/${line}.out.txt | grep SubName`
		DOI=`cat ${FILE}/${line}.out.txt | grep DOI`
		DOMAIN=`cat ${FILE}/${line}.out.txt | grep DOMAIN`
		OS=`cat ${FILE}/${line}.out.txt | grep OS`
		OC=`cat ${FILE}/${line}.out.txt | grep OC`
		InterPro=`cat ${FILE}/${line}.out.txt | grep InterPro`
		TIGRFAMs=`cat ${FILE}/${line}.out.txt | grep TIGRFAMs`

		echo '['`date "+%H:%M:%S"`'] - Retrieve UNIPROT meta-data ';
		echo -e $line '\t' $SubName '\t' $NCBI '\t' $OS '\t' $OC '\t' $EMBL '\t' $KEGG '\t' $PATRIC '\t' $KO '\t' $PFAM  '\t' $DOMAIN '\t' $InterPro  '\t' $TIGRFAMs >> test_uniprot2NCBI.txt


		NCBI2=`echo ${NCBI} | awk -F ';' '{print $2}'`
		
		#echo $NCBI2
		
		if [ -n "$NCBI2" ]; then
    		echo '['`date "+%H:%M:%S"`'] - Download AA sequence ';
			esearch -db protein -query "$NCBI2" < /dev/null | efetch -format fasta >> testFASTA.fasta
			#echo '['`date "+%H:%M:%S"`'] - Retrieve NCBI protein meta-data';
			#esearch -db protein -query "$NCBI2"  | efetch -format xml | xtract -Pattern Seq-entry -first Textseq-id_accession Prot-ref_name_E NCBIeaa >> uni2NCBI.txt 
			#echo '['`date "+%H:%M:%S"`'] - Retrieve sequence specific NCBI taxonomy data';
			#esearch -db protein -query "$NCBI2" | elink -target taxonomy | efetch -format xml | xtract -Pattern Taxon -first TaxId ScientificName ParentTaxId Lineage CreateDate >> tax.UNI2NBCI.txt 
		else
    		echo '['`date "+%H:%M:%S"`'] - NO NCBI-ID LINK FOUND ';
		fi		

	done
	
	
echo '===================================================';


#esearch -db protein -query "WP_038636709.1" | efetch -format fasta
#esearch -db protein -query "WP_038636709.1" | efetch -format xml | xtract -Pattern Seq-entry -first Textseq-id_accession Prot-ref_name_E NCBIeaa >> uni2NCBI.txt 
#esearch -db protein -query "WP_038636709.1" | elink -target taxonomy | efetch -format xml | xtract -Pattern Taxon -first TaxId ScientificName ParentTaxId Lineage CreateDate
#!/bin/sh


FILENAME=$1
echo $FILENAME

FILE=`basename ${FILENAME} .txt`
mkdir ${FILE};

cat ${FILENAME} | while read line
	do
		echo ${line}
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

		echo -e $line '\t' $SubName '\t' $NCBI '\t' $OS '\t' $OC '\t' $EMBL '\t' $KEGG '\t' $PATRIC '\t' $KO '\t' $PFAM  '\t' $DOMAIN '\t' $InterPro  '\t' $TIGRFAMs >> test_uniprot2NCBI.txt
	done
	
	
#cat /Users/sa01fd/testUltra/test_uniprot2NCBI.txt | awk -F '\t' '{print $3}' | awk -F ';' '{print $2}'


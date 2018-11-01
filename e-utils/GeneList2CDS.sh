#from PROTEIN ID TO NCBI
#efetch -db protein -format fasta_cds_na -id XP_003399879.1




#a one liner, if ever needed
#esearch -db protein -query MARHY0005 | elink -target nuccore|efetch -format fasta_cds_na | sed -n -e '/MARHY0005/,/>/ p' | sed '$d'



#download the MARHY genome
#esearch -db protein -query MARHY0005 | elink -target nuccore|efetch -format fasta_cds_na > MARHY.genome





filename='geneList.txt'

cat $filename | while read line
	do
   		cat MARHY.genome | sed -n -e "/$line/,/>/ p" | sed '$d'	>> ExtractedFasta.fasta
	done


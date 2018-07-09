#!/bin/sh
# Protein-Protein BLAST 2.7.1+
# this scripts needs BLAST+ < 2.6


#	MAKEBLASTDB:		
# take all the fasta files in a certain directory and make a combined file to be converted to blastdb

#cat *.fasta > all.fasta

#/Users/sa01fd/Genomics/TEST_enrichM/genome_proteins/Marinobacter_hydrocarbonoclasticus_VT8.faa 

#cat /Users/sa01fd/Genomics/TEST_enrichM/genome_proteins/*.faa > all.fasta
#makeblastdb -in all.fasta  -out protDB_marinobacter -dbtype prot

#makeblastdb -in /Users/frederik/Downloads/GCF_000284615.1_ASM28461v1_protein.faa  -out testDB -dbtype prot

echo '===================================================';
echo 'blastPultra.sh version beta 0.1';
echo '===================================================';
starttime=`date +%s`



echo '['`date "+%H:%M:%S"`'] - Starting blastPultra.sh';



runtime=$((end-start))

#	BLASTP:	
#	LOCAL:	
#If those files are still hiden in a prokka folder ... execute the below
echo '['`date "+%H:%M:%S"`'] - Extracting genome annotations (.faa) from prokka directory';
find prokka -name \*.faa -exec cp {} ./faa \;

#make a gaint fasta-file by merging annotation files (from ./faa/)
echo '['`date "+%H:%M:%S"`'] - Merging annotations for blastdb generation';
cat ./faa/*.faa > all.fasta

#make blast protein database using all.fasta and save as protDB_marinobacter
echo '['`date "+%H:%M:%S"`'] - Making blastdb ...';
makeblastdb -in all.fasta  -out protDB_marinobacter -dbtype prot

#local blast, output saved as results.out
echo '['`date "+%H:%M:%S"`'] - Running blastp ...';
blastp -query ./H8WCU7.fasta -db protDB_marinobacter -out results.out -outfmt "6 qseqid sseqid slen qstart qend length mismatch gapopen gaps sseq" -evalue 1e-20

#Generate a fasta file from it
echo '['`date "+%H:%M:%S"`'] - Saving fasta file ...';
 awk 'BEGIN { OFS = "\n" } { print ">"$2, $10 }' results.out > local_blast.fasta
 
#Generate a file containing all the accession numbers, used to trace back the organism
 awk 'BEGIN { OFS = "\n" } { print $2}' results.out > lookup.txt


#remove files from 'potential' previous run
[ -e Origin.txt ] && rm Origin.txt
[ -e abundance_gene.txt ] && rm abundance_gene.txt

echo '['`date "+%H:%M:%S"`'] - Calculating abundance table ...';
#trace back the organism
cat lookup.txt | while read line
do
   find faa/*.faa -type f -print0 | xargs -0 grep -l $line | xargs -I"{}" basename {} .faa >> Origin.txt
done

#determin abundance of genes, per organisms
cat Origin.txt | while read line
do
	countl=`grep -c ${line} Origin.txt`
	echo  ${line} ${countl} >> abundance_gene.txt
done

#remove duplicate to obtain a clean abundance table
sort -u -k 1,1 abundance_gene.txt > abundance_gene2.txt


paste Origin.txt results.out | column -s $'\t' -t > blast_local_results.tmp
awk 'BEGIN {print "sOrganism\tqseqid\tsseqid\tslen\tqstart\tqend\tlength\tmismatch\tgapopen\tgaps\tsseq"} {print}' blast_local_results.tmp > blast_local_results.csv
rm blast_local_results.tmp

echo '['`date "+%H:%M:%S"`'] - Formatting output table [blast_local_results.csv]...';




endtime=`date +%s`
runtime=$((endtime-starttime))
echo '['`date "+%H:%M:%S"`'] - Finished .... [total runtime:' ${runtime} 'seconds]';
echo '===================================================';





 
 
 
#blastp -query ./aupB/aupB_MARHY.fasta -db testDB -out results.out -outfmt 6 -evalue 0.005
#blastp -query ./aupB/aupB_MARHY.fasta -db testDB -out results.out -outfmt "6 qseqid sseqid slen qstart qend length mismatch gapopen gaps sseq" -evalue 0.005

#	REMOTE:	
#blastp -query ./aupB/aupB_MARHY.fasta -db nr -out results.out -outfmt "6 qseqid sseqid sseq sacc staxids sscinames slen evalue bitscore score qstart qend length mismatch gapopen gaps" -evalue 1e-12 -remote

#blastp -query ./aupB/aupB_MARHY.fasta -db nr -out results.out -outfmt "6 qseqid sseqid sseq sacc staxids sscinames slen evalue bitscore score qstart qend length mismatch gapopen gaps" -evalue 1e-12 -remote -ungapped



#	TO FASTA FILE:	
#awk 'BEGIN { OFS = "\n" } { print ">"$2, $3 }' results.out



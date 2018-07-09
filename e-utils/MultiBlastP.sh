#!/bin/sh
# Protein-Protein BLAST 2.7.1+
# this scripts needs BLAST+ < 2.6


echo '===================================================';
echo 'multiBlastp.sh version beta 0.1';
echo '===================================================';

starttime=`date +%s`
echo '['`date "+%H:%M:%S"`'] - Starting multiBlastp.sh';

MULTIFASTA='salicylate.fasta';
outdir=${MULTIFASTA};
mkdir 'salicylate';

echo '['`date "+%H:%M:%S"`'] - Extracting genome annotations (.faa) from prokka directory';
find prokka -name \*.faa -exec cp {} ./faa \;

#make a gaint fasta-file by merging annotation files (from ./faa/)
echo '['`date "+%H:%M:%S"`'] - Merging annotations for blastdb generation';
cat ./faa/*.faa > all.fasta

#make blast protein database using all.fasta and save as protDB_marinobacter
echo '['`date "+%H:%M:%S"`'] - Making blastdb ...';
makeblastdb -in all.fasta  -out protDB_marinobacter -dbtype prot


echo '...';
echo '...';
echo '...';


echo '['`date "+%H:%M:%S"`'] - Splitting Multifasta file';
awk -F "|" '/^>/ {F = "./salicylate/"$2".fasta"} {print > F}' ${MULTIFASTA};



for f in ./salicylate/*.fasta; do 
	FILENAME=`basename ${f} .fasta`
	echo '-------------------------------------------';
	echo "Processing $FILENAME"
	echo '-------------------------------------------';
	
	
	#local blast, output saved as results.out
	echo '['`date "+%H:%M:%S"`'] - Running blastp ...';
	blastp -query $f -db protDB_marinobacter -out ${FILENAME}.results.out -outfmt "6 qseqid sseqid slen qstart qend length mismatch gapopen gaps sseq" -evalue 1e-20

	echo '['`date "+%H:%M:%S"`'] - Saving fasta file ...';
 	awk 'BEGIN { OFS = "\n" } { print ">"$2, $10 }' ${FILENAME}.results.out > ${FILENAME}.local_blast.fasta
 
	#Generate a file containing all the accession numbers, used to trace back the organism
 	awk 'BEGIN { OFS = "\n" } { print $2}' ${FILENAME}.results.out > ${FILENAME}.lookup.txt

	echo '['`date "+%H:%M:%S"`'] - Calculating abundance table ...';
	#trace back the organism
	cat ${FILENAME}.lookup.txt | while read line
	do
   		find faa/*.faa -type f -print0 | xargs -0 grep -l $line | xargs -I"{}" basename {} .faa >> ${FILENAME}.Origin.txt
	done

	#determin abundance of genes, per organisms
	cat ${FILENAME}.Origin.txt | while read line
	do
		countl=`grep -c ${line} ${FILENAME}.Origin.txt`
		echo ${FILENAME} ${line} ${countl} >> ${FILENAME}.abundance_gene.txt
	done

	#remove duplicate to obtain a clean abundance table
	sort -u -k 1,1 ${FILENAME}.abundance_gene.txt > ${FILENAME}.abundance_gene2.txt

	echo '['`date "+%H:%M:%S"`'] - Formatting output table [blast_local_results.csv]...';

	paste ${FILENAME}.Origin.txt ${FILENAME}.results.out | column -s $'\t' -t > ${FILENAME}.blast_local_results.tmp
	awk 'BEGIN {print "sOrganism\tqseqid\tsseqid\tslen\tqstart\tqend\tlength\tmismatch\tgapopen\tgaps\tsseq"} {print}' ${FILENAME}.blast_local_results.tmp > ${FILENAME}.blast_local_results.csv
	rm ${FILENAME}.blast_local_results.tmp
done
echo '===================================================';


for gene2 in ./*abundance_gene2.txt; do
		echo "Processing $gene2"
		cat $gene2 >> Long_abundance.txt
done


echo '['`date "+%H:%M:%S"`'] - Formatting overall abundance table...';
Rscript -e 'loadtab<-read.table("~/testUltra/Long_abundance.txt");wide = reshape(loadtab, idvar = "V1", timevar = "V2", direction = "wide");names(wide) <- gsub("V3.", "", names(wide));rownames(wide)=wide$V1;wide=wide[,2:ncol(wide)];write.table(wide,file="abundance_table.txt",quote=FALSE)'


#read.table('~/testUltra/abundance_table.txt')


endtime=`date +%s`
runtime=$((endtime-starttime))
echo '['`date "+%H:%M:%S"`'] - Finished .... [total runtime:' ${runtime} 'seconds]';
echo '===================================================';










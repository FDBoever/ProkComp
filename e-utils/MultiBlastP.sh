#!/bin/sh
# Protein-Protein BLAST 2.7.1+
# this scripts needs BLAST+ < 2.6

usage()
{
    echo "usage: sysinfo_page [[[-f --fasta ] [-e --evalue]] | [-h --help]]"
}

#SET DEFAULTS
#interactive=
evalue=12

while [ "$1" != "" ]; do
    case $1 in
        -f | --fasta )           shift
                                fasta=$1
                                ;;
        #-d | --blastdb )    	shift
        #						blastdb=$1
        #                        ;;
        -e | --evalue )         shift
        						evalue=$1
                                ;;
		-h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

echo '===================================================';
echo 'multiBlastp.sh version beta 0.1';
echo '===================================================';

starttime=`date +%s`
echo '['`date "+%H:%M:%S"`'] - Starting multiBlastp.sh';

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


echo ''
echo ' ::: INPUT :::'

echo 'FastaFile:' ${fasta}
echo 'e-value:' ${evalue}

FILENAME=`basename ${fasta} .fasta`
outdir=${FILENAME};
mkdir ${FILENAME};
mkdir ${FILENAME}'/InputFasta/';
mkdir ${FILENAME}'/SplitFasta/';

SPLITPATH=${FILENAME}"/SplitFasta/"
echo $SPLITPATH
cp ${fasta} ${FILENAME}'/InputFasta/'

echo '['`date "+%H:%M:%S"`'] - Splitting Multifasta file';
#awk -F "|" '/^>/ {F = "./"$SPLITPATH"/"$2".fasta"} {print > F}' ${fasta};
awk -v var="${SPLITPATH}/" -F "|" '/^>/ {F = var $2".fasta"} {print > F}' ${fasta};
#"./"${FILENAME}"/"




for fas in $SPLITPATH*.fasta; do 
	fname=`basename ${fas} .fasta`
	echo '-------------------------------------------';
	echo "Processing $fname"
	echo '-------------------------------------------';
	
	
	#local blast, output saved as results.out
	echo '['`date "+%H:%M:%S"`'] - Running blastp ...';
	blastp -query $fas -db protDB_marinobacter -out ${FILENAME}/${fname}.results.out -outfmt "6 qseqid sseqid slen qstart qend length mismatch gapopen gaps sseq" -evalue 1e-30

	echo '['`date "+%H:%M:%S"`'] - Saving fasta file ...';
 	awk 'BEGIN { OFS = "\n" } { print ">"$2, $10 }' ${FILENAME}/${fname}.results.out > ${FILENAME}/${fname}.local_blast.fasta
 
	#Generate a file containing all the accession numbers, used to trace back the organism
 	awk 'BEGIN { OFS = "\n" } { print $2}' ${FILENAME}/${fname}.results.out > ${FILENAME}/${fname}.lookup.txt

	echo '['`date "+%H:%M:%S"`'] - Calculating abundance table ...';
	#trace back the organism
	cat ${FILENAME}/${fname}.lookup.txt | while read line
	do
   		find faa/*.faa -type f -print0 | xargs -0 grep -l $line | xargs -I"{}" basename {} .faa >> ${FILENAME}/${fname}.Origin.txt
	done

	#determin abundance of genes, per organisms
	cat ${FILENAME}/${fname}.Origin.txt | while read line
	do
		countl=`grep -c ${line} ${FILENAME}/${fname}.Origin.txt`
		echo ${fname} ${line} ${countl} >> ${FILENAME}/${fname}.abundance_gene.txt
	done

	#remove duplicate to obtain a clean abundance table
	sort -u ${FILENAME}/${fname}.abundance_gene.txt > ${FILENAME}/${fname}.abundance_gene2.txt

	echo '['`date "+%H:%M:%S"`'] - Formatting output table [blast_local_results.csv]...';

	paste ${FILENAME}/${fname}.Origin.txt ${FILENAME}/${fname}.results.out | column -s $'\t' -t > ${FILENAME}/${fname}.blast_local_results.tmp
	awk 'BEGIN {print "sOrganism\tqseqid\tsseqid\tslen\tqstart\tqend\tlength\tmismatch\tgapopen\tgaps\tsseq"} {print}' ${FILENAME}/${fname}.blast_local_results.tmp > ${FILENAME}.blast_local_results.csv
	rm ${FILENAME}/${fname}.blast_local_results.tmp
done
echo '===================================================';


for gene2 in ${FILENAME}/*abundance_gene2.txt; do
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










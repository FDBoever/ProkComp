
#!/bin/bash

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
echo 'remoteMultiBlastp.sh version beta 0.1';
echo '===================================================';
echo '['`date "+%H:%M:%S"`'] - Starting multiBlastp.sh';
starttime=`date +%s`
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
	#echo '	blastp -query '$fas' -db nr -out '${FILENAME}'/'${fname}'.results.out -outfmt "6 qseqid sseqid sseq sacc staxids sscinames slen evalue bitscore score qstart qend length mismatch gapopen gaps" -evalue 1e-12 -remote
#'
	#blastp -query $fas -db nr -out ${FILENAME}/${fname}.results.out -outfmt "6 qseqid sseqid sseq sacc staxids sscinames slen evalue bitscore score qstart qend length mismatch gapopen gaps" -evalue 1e-12 -remote

	
	#FIX DA ET IN DE JUSTE MAP KOMT FFS
	
	awk 'BEGIN { OFS = "\n" } { print $2}' ${FILENAME}/${fname}.results.out >>  ${FILENAME}/${fname}.accesion.txt

	taxIDlong=`awk -F"|" '{print $2}' ${FILENAME}/${fname}.accesion.txt | tr '\n' ','`\r
	
	echo "${taxIDlong%?}"
	efetch -db taxonomy -id ${taxIDlong}  -format xml | xtract -pattern Taxon -first ScientificName >> names.txt

	
done
echo '===================================================';


#	REMOTE:	
#blastp -query ./aupB/aupB_MARHY.fasta -db nr -out results.out -outfmt "6 qseqid sseqid sseq sacc staxids sscinames slen evalue bitscore score qstart qend length mismatch gapopen gaps" -evalue 1e-12 -remote
#blastp -query ./aupB/aupB_MARHY.fasta -db nr -out results.out -outfmt "6 qseqid sseqid sseq sacc staxids sscinames slen evalue bitscore score qstart qend length mismatch gapopen gaps" -evalue 1e-12 -remote -ungapped

#	TO FASTA FILE:	
#awk 'BEGIN { OFS = "\n" } { print ">"$2, $3 }' results.out


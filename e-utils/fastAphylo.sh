#!/bin/sh

#	
#  	Java SE Development Kit 10.0.1
#


#	PASTE:		Combine into one; paste several fasta files underneath eachother
#cat *.fasta >> aupA_combined.fasta

usage()
{
    echo "usage: sysinfo_page [[[-f --fasta ]"
}

#SET DEFAULTS
#interactive=
evalue=12

while [ "$1" != "" ]; do
    case $1 in
        -f | --fasta )           shift
                                fasta=$1
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
echo 'fasta2phylo.sh version beta 0.1';
echo '===================================================';

starttime=`date +%s`
echo '['`date "+%H:%M:%S"`'] - Starting fasta2phylo.sh';

echo ''
echo ' ::: INPUT :::'

echo 'FastaFile:' ${fasta}

FILENAME=`basename ${fasta} .fasta`
outdir=${FILENAME};
mkdir ${FILENAME};

echo '['`date "+%H:%M:%S"`'] - Removing unwanted parts in the header: eg UNIPROT:';
sed -i -e 's/UNIPROT\://g' ${fasta} 

#	MUSCLE:	Run default MUSCLE ALIGNMENT using muscle 3.8.31
echo '['`date "+%H:%M:%S"`'] - ALIGN with MUSCLE 3.8.31:';
muscle -in ${fasta} -out ${FILENAME}.aligned.fasta

#	BMGE:		
echo '['`date "+%H:%M:%S"`'] - TRIM alignment positions using BMGE 1.12:';
java -jar /Users/sa01fd/Genomics/BMGE-1.12/BMGE.jar -t AA -g 0.3 -i ${FILENAME}.aligned.fasta -of ${FILENAME}.aligned.fasta.bmge

#	FastTree:	Infer phylogeny using FastTree
echo '['`date "+%H:%M:%S"`'] - Infer phylogeny using FastTre:';
FastTree ${FILENAME}.aligned.fasta.bmge > FastTree_${FILENAME}


endtime=`date +%s`
runtime=$((endtime-starttime))
echo '['`date "+%H:%M:%S"`'] - Finished .... [total runtime:' ${runtime} 'seconds]';
echo '===================================================';
echo ''








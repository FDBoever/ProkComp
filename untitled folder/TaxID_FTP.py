#!/usr/bin/python
import urllib2
import os
import sys
import time
import subprocess


if len(sys.argv) != 3:
    print "USAGE: TaxID2Tax.py <Taxon_id_list> <output path>"
    sys.exit(1)

inputFile = open(sys.argv[1], "rt")
content = inputFile.read()
inputFile.close

num_lines = len(open(sys.argv[1]).readlines(  ))
print(num_lines)



content=''
counter=0
for id in open(sys.argv[1]):
    counter += 1
    if len(id) > int(6):
        #print(str(counter)+"/"+str(num_lines)+" // Fetching %s..." % str(id))
        content = content + '\n' + id
    #else:
    #    print(str(counter)+"/"+str(num_lines)+" ------- FAILED %s..." % str(id))

content = content.replace("\n", ",")


os.system('mkdir ' + sys.argv[2])
os.system('efetch -db taxonomy -id ' + content  + '  -format xml | xtract -pattern Taxon -first ScientificName>>' + sys.argv[2]+'/Taxa.txt')


num_lines2 = len(open(sys.argv[2]+'/Taxa.txt').readlines(  ))
counter=0
for taxon in open(sys.argv[2]+'/Taxa.txt'):
    taxon = taxon.strip()
    counter += 1
    if len(taxon) > int(6):
        os.system('esearch -db taxonomy -query "' + taxon + '[Organism]" | elink -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq>>'+ sys.argv[2]+'/out_ftp.txt')
        
        #test = os.system('esearch -db taxonomy -query "' + taxon + '[Organism]" | elink -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq')

        
    	
        
        
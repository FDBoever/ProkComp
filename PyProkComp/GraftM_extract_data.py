import os
import numpy as np
import shutil
import pandas as pd




inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/GraftM_output/'
outputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/GraftM_extract/'

if not os.path.exists(outputdir):
    os.makedirs(outputdir)

print("")
print("")
print("")
print("")

print('======================================================================================')
print("      Welcome to " + "\033[91m" + "GraftM_extract.py v0.3 \033[0m - latest modifications (02-Aug-2017)")
print("       This script will read GraftM output files from directory <inputdir>")
print("     extract the aligned sequences from all files and compile a new fasta file")
print("  Make sure that the numpy and pandas packages (Python) are installed on your machine")
print("  ")
print("  ")
print("  ")
print("                            contact \033[92m Frederik.deboever@sams.ac.uk \033[0m for more information")
print('======================================================================================')
print(" ")
print(" ")

n_files = 0
for selectedfile in os.listdir(inputdir):
   if selectedfile != ".DS_Store":
       filename = str(inputdir) + str(selectedfile)
       n_files += 1
       print(filename)
print(" ")

print('--------------------------------------------')
print( "           \033[91m" + str(n_files) + "\033[0m marker genes are loaded")
print('-------------------------------------------')


n_genomes = 0
for selectedfile in os.listdir(inputdir):
   if selectedfile != ".DS_Store":
       filename = str(inputdir) + str(selectedfile)
       n_genomes = 0
       
       fastaString = ""
        
       for selectedGenome in os.listdir(filename):
           if os.path.isdir(str(inputdir) + str(selectedfile) + "/" + str(selectedGenome)):
               n_genomes += 1
               #print(str(inputdir) + str(selectedfile) + "/" + str(selectedGenome))
               
               
               for file in os.listdir(str(inputdir) + str(selectedfile) + "/" + str(selectedGenome) + "/"):
                   if "_hits.aln.fa" in file:
                                            
                       
                       with open(str(inputdir) + str(selectedfile) + "/" + str(selectedGenome) + "/" + file) as f:
    		               content = f.readlines()
    		               content = [x.strip() for x in content]
    		               fastaString += "\n>" + str(selectedGenome) + "\n" + content[1]

    		               if len(content) ==4:
    		                   fastaString += "\n>" + str(selectedGenome) + "\n" + content[3]
    		               if len(content) ==6:
    		                   fastaString += "\n>" + str(selectedGenome) + "\n" + content[5]

        
       print(fastaString)
       text_file = open(outputdir + selectedfile + ".fasta", "w")
       text_file.write(fastaString)
       text_file.close()              
               
               
       print(" ")
       print('--------------------------------------------')
       print( "    \033[91m" + str(n_genomes) + "\033[0m detected genomes for " + str(selectedfile))
       print('-------------------------------------------')
       print(" ")        
   
   
   
       
       
print(" ")


print(' ')
print('======================================================================================')
print('\033[92m saving output file(s)\033[0m')
print(">>>  /Users/sa01fd/DATA/MarinobacterGenomics/GraftM/GraftM_extract/")
print(' ')
print("\033[91m DONE!!! \033[0m   Thanks for using GraftM_extract()")
print('======================================================================================')
print(' ')

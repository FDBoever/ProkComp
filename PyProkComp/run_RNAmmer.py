import os
import shutil


inputdir = '/Users/sa01fd/DATA/MarinobacterGenomics/FNA_Marinobacter_restricted/'
outputstring = '/Users/sa01fd/DATA/MarinobacterGenomics/RNANMMER/' 

string = ''



for selectedfile in os.listdir(inputdir):
   filename = str(inputdir) + str(selectedfile)
   outname = str(outputstring) + str(selectedfile)
   print(filename)
   print('-------------------------------------------')
   print(' Proccess: ' + filename) 
   print('-------------------------------------------')
   string += 'rnammer -S bac -h ' + selectedfile + '.hmm_report.txt -f ' + selectedfile + '.fasta ' + filename + '; '
print(string)







 

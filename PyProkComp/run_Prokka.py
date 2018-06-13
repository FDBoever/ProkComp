import os
import shutil
inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/Assemblies/Marinobacter2/'

string = ''

for selectedfile in os.listdir(inputdir):
   filename = str(inputdir) + str(selectedfile)
   #outname = str(outputstring) + str(selectedfile)
   print(filename)
   print('-------------------------------------------')
   print(' Proccess: ' + filename) 
   print('-------------------------------------------')
   string += 'prokka --kingdom Bacteria --outdir ./PROKKA_output/Prokka_' + str(selectedfile) + ' --genus Marinobacter --locustag ' + str(selectedfile) + ' '  + filename + ' ; ' 
   BashCommand = 'prokka --kingdom Bacteria --outdir ./PROKKA_output/Prokka_' + str(selectedfile) + ' --genus Marinobacter --locustag ' + str(selectedfile) + ' '  + filename + ' ; '
   os.system(BashCommand)
print(string)

import os
import shutil


inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/Assemblies/'
outdir = '/media/data/SANDBOXES/sa01fd/ProkComp/BUSCO_output/' 

string = ''
if not os.path.exists(outdir):
    os.makedirs(outdir)

for selectedfile in os.listdir(inputdir):
   filename = str(inputdir) + str(selectedfile)
   #outname = str(outputstring) + str(selectedfile)
   outputPath = "/media/data/SANDBOXES/sa01fd/ProkComp/BUSCO_output/" + str(selectedfile[:-4])
   print(filename)
   print('-------------------------------------------')
   print(' Proccess: ' + filename) 
   print('-------------------------------------------')
   if not os.path.exists(outputPath):
       os.makedirs(outputPath)
   string += 'python /opt/softwares/BUSCO/BUSCO_v1.1b1.py -in ' + filename + ' -o ' + outputPath + ' -l bacteria -m genome; '
   
   
print(string)



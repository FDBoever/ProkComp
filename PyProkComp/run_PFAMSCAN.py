import os
import shutil


inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/Annotations/Marinobacter/'
outdir = '/media/data/SANDBOXES/sa01fd/ProkComp/PfamScan_output/Marinobacter/' 

string = ''



#check whether the output directory (Annotations) exists and make accordingly
#filenames are trimmed, getting rid of the run date and ".fna"

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
    
for selectedfile in os.listdir(inputdir):
   filename = str(inputdir) + str(selectedfile)
   outname = str(outdir) + str(selectedfile)
   print(filename)
   print('-------------------------------------------')
   print(' Proccess: ' + filename) 
   print('-------------------------------------------')
   BashCommand = './PfamScan/pfam_scan.pl -fasta ' + filename + ' -outfile ' + outname + '.txt -dir ./PfamScan/ '  + '-cpu 4; '

   print(BashCommand)
   string += './PfamScan/pfam_scan.pl -fasta ' + filename + ' -outfile ' + outname + '.txt -dir ./PfamScan/ '  + '-cpu 4; '
   os.system(BashCommand)

print(string)







 

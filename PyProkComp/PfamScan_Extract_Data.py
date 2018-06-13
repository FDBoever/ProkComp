import os
import numpy as np
import shutil
import pandas as pd


inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/PfamScan_output/Marinobacter/'
outdir = '/media/data/SANDBOXES/sa01fd/ProkComp/PfamScan_output/'
print("")
print("")
print("")
print("")

print('======================================================================================')
print("      Welcome to " + "\033[91m" + "PfamScan_extract.py v0.3 \033[0m - latest modifications (22-08-2017)")
print("       This script will read PfamScan output files from directory <inputdir>")
print("     extract the pfam domain columns from all files and compile a new dataframe")
print("  Make sure that the numpy and pandas packages (Python) are installed on your machine")
print("  ")
print("  \033[95mPfamScan\033[0m (Sanger tools) searches FASTA sequences against a library of pfam HMM ")
print("  This tool can be obtained from EMBL website; and offline installation is recommended")
print("      a paper to cite is: \033[95m Mistry et al. 2007; BMC Bioinformatics. 2007; 8; 298 \033[0m")
print("  ")
print("  ")
print("                            contact \033[92m Frederik.deboever@sams.ac.uk \033[0m for more information")
print('======================================================================================')
print(" ")
print(" ")

n_files = 0
for selectedfile in os.listdir(inputdir):
   filename = str(inputdir) + str(selectedfile)
   num_lines = sum(1 for line in open(filename))
   n_files += 1
   print("opening >> \033[92m " + str(num_lines) + "  \033[0m lines : \033[92m" + selectedfile + "\033[0m")
print(" ")

print('--------------------------------------------')
print( "           \033[91m" + str(n_files) + "\033[0m files are loaded")
print('-------------------------------------------')

#df = np.array(["seq id","alignment start","alignment end","envelope start","envelope end","hmm acc","hmm name","type","hmm start","hmm end","hmm length","bit score","E-value","significance","clan"])

df = pd.DataFrame(columns = ("genome","domain","Clan"))
print(df)

#df = np.array(['Marinobacter_zhejiangensis_CGMCC1_7061.gbk.fna_03504', '571', '655', '570', '657', 'PF00679.22', 'EFG_C', 'Domain', '2', '87', '89', '104.6', '1.9e-30', '1', 'CL0437'])
print(" ")
total_index = 0

n_files = 0
for selectedfile in os.listdir(inputdir):

	filename = str(inputdir) + str(selectedfile)
	n_files += 1
	print("analysing file  [" + str(n_files) + "] >>> \033[92m" + selectedfile + "\033[0m")
   
	with open(filename, "r+") as f:
		n_ommited = 0
		n_lines = 0
		for line in f:
			if "#" in line:
				#print(line)
				n_ommited += 1
			elif ".fna" in line:
				#print(line)
				anl_line = line
				df.loc[total_index] = [selectedfile,anl_line.split()[6],anl_line.split()[14]]
				#print(anl_line.split()[6])
				n_lines += 1
				total_index += 1
		print("\033[92m" + str(n_ommited)+ "\033[0m lines were ommitted line due to masking (#)")
		print("\033[92m" + str(n_lines-1)+ "\033[0m pfam domains (PfamScan)")
print(" ")
print(" ")
print('--------------------------------------------')
print("Data from multiple files were extracted and fit into a new matrix")
print(" ")
print(df)


df.to_csv( outdir + "PfamScan_output_oceano.csv", sep='\t')
print(' ')
print('======================================================================================')
print('\033[92m saving output file(s)\033[0m')
print(">>>  " + outdir + "PfamScan_output.csv")
print(' ')
print("\033[91m DONE!!! \033[0m   Thanks for using PfamScan_extract()")
print('======================================================================================')
print(' ')

import os
import numpy as np
import shutil
import pandas as pd


inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/'
outputdir = ''
print("")
print("")
print("")
print("")
print('======================================================================================')
print("      Welcome to " + "\033[91m" + "CheckM_filter.py \033[0m - latest modifications (21-08-2017)")
print("       This script will read CheckM outputfile  <inputfile>")
print("     and copy all the assemblies of > 98% completeness to a new folder QC_Assemblies")
print("  Make sure that the numpy and pandas packages (Python) are installed on your machine")
print("  ")
print("                            contact \033[92m Frederik.deboever@sams.ac.uk \033[0m for more information")
print('======================================================================================')
print(" ")
print(" ")


#Reads CheckM output file and obtains some file statistics
filename = str(inputdir) + "/CheckM_output/qa2.txt"
num_lines = sum(1 for line in open(filename))
print("opening >> \033[92m " + str(num_lines) + "  \033[0m lines : \033[92m" + filename + "\033[0m")
print(" ")

#setting variables that correspond to panda matrices to be filled in script
df = pd.DataFrame(columns = ("genome","completeness","QC"))
dfommited = pd.DataFrame(columns = ("genome","completeness"))
dfkeep = pd.DataFrame(columns = ("genome","completeness"))

#Assigning basic count variables
total_index = 0
n_files = 0

#the threshold for genome completeness 
threshold = 98.00

if not os.path.exists(inputdir + "QC_Assemblies/"):
    os.makedirs(inputdir + "QC_Assemblies/")

#For each line in the file, check criteria, and copy the respective information in the matrices
with open(filename, "r+") as f:
		n_ommited = 0
		n_lines = 0
		for line in f:
			if "---" in line:
				#print(line)
				n_ommited += 1
			elif "Bin" in line:
				n_ommited += 1
			else:
				#print(line)
				anl_line = line
				if float(anl_line.split()[6]) >= threshold:
					QC = "yes"
					dfkeep.loc[total_index] = [anl_line.split()[0],anl_line.split()[6]]
					print( "Copying genome: \033[92m" + anl_line.split()[0] + "\033[91m (" + anl_line.split()[6] + " %) \033[0m")
					shutil.copy2("./Assemblies/" + anl_line.split()[0] + ".fna", "./QC_Assemblies/" + anl_line.split()[0] + ".fna")
				else:
					QC = "no"
					dfommited.loc[total_index] = [anl_line.split()[0],anl_line.split()[6]]
				df.loc[total_index] = [anl_line.split()[0],anl_line.split()[12],QC]
				#print(anl_line.split()[6])
				n_lines += 1
				total_index += 1
#print(" ")
#print(" ")
#print('--------------------------------------------')
#print("Data from multiple files were extracted and fit into a new matrix")
#print(" ")

print('--------------------------------------------')
print("\033[91m FOLLOWING GENOMES ARE OMMITED IN FURTHER ANALYSIS\033[0m")
print(dfommited)
print('--------------------------------------------')

print('--------------------------------------------')
print("\033[92m FOLLOWING GENOMES ARE KEPT FOR FURTHER ANALYSIS\033[0m")
print(dfkeep)
print('--------------------------------------------')



#df.to_csv("/Users/sa01fd/Downloads/PfamScan_output.csv", sep='\t')
print(' ')
print('======================================================================================')
print('\033[92m saving output file(s)\033[0m')
print(">>>  /Users/sa01fd/Genomics/ProkComp/CheckM.txt")
print(' ')
print("\033[91m DONE!!! \033[0m   Thanks for using CheckM_filter)")
print('======================================================================================')
print(' ')

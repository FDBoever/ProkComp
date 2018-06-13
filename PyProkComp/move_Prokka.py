import os
import numpy as np
import shutil
import pandas as pd


inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/Annotations/Marinobacter2/PROKKA_output/'
outdir = '/media/data/SANDBOXES/sa01fd/ProkComp/Annotations/Marinobacter2/'
print("")
print("")
print("")
print("")

print('======================================================================================')
print("      Welcome to " + "\033[91m" + "move_PROKKA.py v0.3 \033[0m - latest modifications (22-Aug-2017)")
print("       This script will copy PROKKA output files from directory <inputdir> to <outdir>")
print("  ")
print("  ")
print("                            contact \033[92m Frederik.deboever@sams.ac.uk \033[0m for more information")
print('======================================================================================')
print(" ")
print(" ")


#check whether the output directory (Annotations) exists and make accordingly
#filenames are trimmed, getting rid of the run date and ".fna"

if not os.path.exists(outdir):
    os.makedirs(outdir)



#look in PROKKA output structure, in each genome folder for files that end with .faa
for basename in os.listdir(inputdir):
    pathname = os.path.join(inputdir, basename)   
    
    if os.path.isdir(pathname):
        for file in os.listdir(pathname):
            if file.endswith('.faa'):
                pathname2 = os.path.join(pathname, file)
                if os.path.isfile(pathname2):
                    shutil.copy2(pathname2, os.path.join(outdir, pathname[85:-4] + ".faa"))
                    print(pathname[85:-4])  
                    print("Annotation copied to \033[96m " + os.path.join(outdir, file[:-17] + ".faa \033[0m"))
       
print(" ")







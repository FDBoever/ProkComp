import os
import shutil
import sys
import numpy as np
import shutil
import pandas as pd

inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/'


################################################################
#
# PIPELINE FOR PROKARYOTIC GENOME COMPARISON
#
#	How to run it:
#		... deposit *.fna in unique directory (name = <AnalysisName>) in ./Assemblies/
#		... Change directory to ./ProcComp/ and run the python script as follows: 
#			... python ./ProkComp.py <AnalysisName>
#
#################################################################

print('======================================================================================')
print("      Welcome to " + "\033[96m" + "ProkComp.py v1 \033[0m - latest modifications (23-Aug-2017)")
print("                            contact \033[95m Frederik.deboever@sams.ac.uk \033[0m for more information")
print('======================================================================================')
print(" ")
print(" ")

#setting environment to recognise Prodigal and PPlacer
os.system("export PATH=$PATH:/opt/softwares/prodigal/2.6.3")
os.system("export PATH=$PATH:/opt/softwares/pplacer/v1.1.alpha17/")
os.system("export PATH=$PATH:/media/data/SANDBOXES/sa01fd/ProkComp/SingalP/")


AnalysisName = sys.argv[1]
os.system("echo Your analysis is called \033[96m" + AnalysisName + "\033[0m")
os.system("echo ... no corrupt files are detected ")
os.system("echo ... data and output files are stored under \033[96m" + AnalysisName + "\033[0m subdirectories")
os.system("echo ... \033[96m Prodigal \033[0m accesible from environment")
os.system("echo ... \033[96m pplacer \033[0m accesible from environment")
os.system("echo ... \033[96m SignalP \033[0m accesible from environment")
print(" ")
print(" ")


################################################################
#
# RUNNING CHECKM
#
#################################################################


print('**********************************************')
print(" Step 1: QUALITY CONTROL  using \033[92m CheckM \033[0m" )
print('**********************************************')


print(" ")
print(" ")
outputdir = inputdir + "/CheckM_output/" + AnalysisName + "/"
outdir = inputdir + "/CheckM_output/"

if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(outputdir):
    os.system("echo ... generating " + outputdir)
    os.makedirs(outputdir)

os.system("echo ... files will be stored in  " + outputdir)
os.system("rm -r ./CheckM_output/" + AnalysisName )

#Running the checkM pipeline using 16 threads
os.system("checkm tree -t 16 " + inputdir + "/Assemblies/" + AnalysisName + " " + inputdir + "/CheckM_output/" + AnalysisName + "/tree/") 
os.system("checkm tree_qa " + inputdir + "/CheckM_output/" + AnalysisName + "/tree/") 
os.system("checkm lineage_set " + inputdir + "/CheckM_output/" + AnalysisName + "/tree/ " + inputdir + "/CheckM_output/" + AnalysisName + "/tree/marker_file.txt") 
os.system("checkm analyze " + inputdir + "/CheckM_output/" + AnalysisName + "/tree/marker_file.txt " + inputdir +"/Assemblies/" + AnalysisName + " " + inputdir + "/CheckM_output/" + AnalysisName + "/") 
os.system("checkm qa -o 2 " + inputdir + "/CheckM_output/" + AnalysisName + "/tree/marker_file.txt " + inputdir +"/CheckM_output/" + AnalysisName + "/ -f " + inputdir + "/CheckM_output/" + AnalysisName + "/qa2.txt")



################################################################
#
# CHECKM FILTERING
#
#################################################################

#Reads CheckM output file and obtains some file statistics
filename = str(inputdir) + "/CheckM_output/" + AnalysisName + "/qa2.txt"
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

if not os.path.exists(inputdir + "QC_Assemblies/"  + AnalysisName + "/"):
    os.makedirs(inputdir + "QC_Assemblies/" + AnalysisName + "/")

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
					shutil.copy2("./Assemblies/" + AnalysisName + "/" + anl_line.split()[0] + ".fna", "./QC_Assemblies/" + AnalysisName + "/" + anl_line.split()[0] + ".fna")
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

################################################################
#
# PROKKA: Genome annotation 
#
#################################################################

inputdir2 = '/media/data/SANDBOXES/sa01fd/ProkComp/QC_Assemblies/' + AnalysisName + "/"


print(" ")
print(" ")
print('**********************************************')
print(" Step 3: GENOME ANNOTATION  using \033[92m PROKKA \033[0m" )
print('**********************************************')
print(" ")
print(" ")

string = ''

for selectedfile in os.listdir(inputdir2):
   filename = str(inputdir2) + str(selectedfile)
   print(filename)
   print('-------------------------------------------')
   print(' Proccess: ' + filename) 
   print('-------------------------------------------')
   string += 'prokka --kingdom Bacteria --outdir ./PROKKA_output/' + AnalysisName + '/Prokka_' + str(selectedfile) + ' --genus Marinobacter --locustag ' + str(selectedfile) + ' '  + filename + ' --gram neg --cpus 16 --force; ' 
   BashCommand = 'prokka --kingdom Bacteria --outdir ./PROKKA_output/' + AnalysisName + '/Prokka_' + str(selectedfile) + ' --genus Marinobacter --locustag ' + str(selectedfile) + ' '  + filename + ' --gram neg --cpus 16 --force; '
   os.system(BashCommand)
print(string)



################################################################
#
# MOVE PROKKA FILES (PROKKA NEEDS TO BE INSTALLED PROPERLY)
#
#################################################################

print(" ")
print(" ")
print('**********************************************')
print(" Step 4: STORING ANNOTATIONS in folder" )
print('**********************************************')
print(" ")
print(" ")

inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/PROKKA_output/' + AnalysisName + "/"
outdir = '/media/data/SANDBOXES/sa01fd/ProkComp/Annotations/' + AnalysisName + "/"

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
                    shutil.copy2(pathname2, os.path.join(outdir, file[:-17] + ".faa"))  
                    print("Annotation copied to \033[96m " + os.path.join(outdir, file[:-17] + ".faa \033[0m"))       
print(" ")


################################################################
#
# BUSCO (NEEDS TO BE INSTALLED PROPERLY)
#
#################################################################

inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/Annotations/' + AnalysisName + '/'
outdir = '/media/data/SANDBOXES/sa01fd/ProkComp/BUSCO_output/' + AnalysisName + '/'

print(" ")
print(" ")
print('**********************************************')
print(" Step 2: QUALITY CONTROL  using \033[92m BUSCO \033[0m" )
print('**********************************************')
print(" ")
print(" ")

string = ''
if not os.path.exists(outdir):
    os.makedirs(outdir)

os.chdir('./BUSCO_output/')

for selectedfile in os.listdir(inputdir):
   filename = str(inputdir) + str(selectedfile)
   outputPath = outdir + str(selectedfile[:-4])
   print(filename)
   print('-------------------------------------------')
   print(' Proccess: ' + filename) 
   print('-------------------------------------------')
   if not os.path.exists(outputPath):
       os.makedirs(outputPath)
   #string += 'python /opt/softwares/BUSCO/BUSCO_v1.1b1.py -in ' + filename + ' -o ' + outputPath + ' -l bacteria -m genome; '
   os.system('python /opt/softwares/BUSCO/v3/scripts/run_BUSCO.py -i ' + filename + ' -o ' + str(selectedfile[:-4]) + ' -l bacteria_odb9 -m prot -f;')
   
os.chdir('..')

outputstring = "genome \t S \t D \t F \t M \n"
count = 0

for DIR in os.listdir('/media/data/SANDBOXES/sa01fd/ProkComp/BUSCO_output/'):
    filename = '/media/data/SANDBOXES/sa01fd/ProkComp/BUSCO_output/' + str(DIR)
    if "run_" + AnalysisName in DIR:
        with open(filename + "/short_summary_" + DIR[4:] + ".txt") as f:
            lines = f.readlines()
            outputstring += DIR[4:] + "\t" + str(lines[10].split('\t')[1]) + "\t"+ str (lines[11].split('\t')[1]) + "\t"+ str (lines[12].split('\t')[1]) + "\t"+ str (lines[13].split('\t')[1]) + "\n"
            #content = f.readlines()
    	    #content = [x.strip() for x in content]
    	    #fastaString += "\n>" + str(selectedGenome) + "\n" + content[1]
            count += 1
print(outputstring)

text_file = open(outdir + AnalysisName + "_BUSCO_SUMMARY.txt", "w")
text_file.write(outputstring)
text_file.close()      


################################################################
#
# GraftM for the detection of ribosomal protein coding genes
#
#################################################################


print(" ")
print(" ")
print('**********************************************')
print(" Step 5: Ribosomal gene detection using \033[92m GraftM \033[0m" )
print('**********************************************')
print(" ")
print(" ")

if not os.path.exists('/media/data/SANDBOXES/sa01fd/ProkComp/GraftM_output/' + AnalysisName + "/"):
    os.makedirs('/media/data/SANDBOXES/sa01fd/ProkComp/GraftM_output/' + AnalysisName + "/")


os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.01.2013_08_greengenes_61_otus.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/test_16S;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.07.ribosomal_protein_L2_rplB.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.07_L2_rplB;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.08.ribosomal_protein_L3_rplC.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.08_L3_rplC;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.09.ribosomal_protein_L5_rplE.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.09_L5_rplE;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.10.ribosomal_protein_L6_rplF.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.10_L6_rplF;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.11.ribosomal_protein_L10.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.11_L10;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.12.ribosomal_protein_L11_rplK.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.12_L11_rplK;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.13.ribosomal_protein_L14b_L23e_rplN.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.13_L14b_L23e_rplN;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.14.ribosomal_protein_L16_L10E_rplP.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.14_L16_L10E_rplP;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.15.ribosomal_protein_S2_rpsB.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.14_S2_rpsB;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.16.ribosomal_protein_S5.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.16_S5;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.17.ribosomal_protein_S7.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.17_S7;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.18.ribosomal_protein_S10_rpsJ.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.18_S10_rpsJ;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.19.ribosomal_protein_S12_S23.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.19_S12_S23;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.20.ribosomal_protein_S15P_S13e.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.20_S15P_S13e;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.21.ribosomal_protein_S19_rpsS.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.21_S19_rpsS;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.27.methyl_coenzyme_reductase_alpha_subunit.mcrA.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.27_methyl_coenzyme_reductase_alpha_subunit_mcrA;')
os.system('graftM graft --forward ./QC_Assemblies/'+ AnalysisName + '/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.40.2013_08_greengenes_97_otus.with_euks.gpkg/ --output_directory ./GraftM_output/'+ AnalysisName + '/o_4.40.2013.8.greengenes;')


inputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/GraftM_output/' + AnalysisName + '/'
outputdir = '/media/data/SANDBOXES/sa01fd/ProkComp/GraftM_extract/' + AnalysisName + '/'

if not os.path.exists(outputdir):
    os.makedirs(outputdir)

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
        
       text_file = open(outputdir + selectedfile + ".fasta", "w")
       text_file.write(fastaString)
       text_file.close()              
                  
       print(" ")
       print('--------------------------------------------')
       print( "    \033[91m" + str(n_genomes) + "\033[0m detected genomes for " + str(selectedfile))
       print('-------------------------------------------')
       print(" ")        
       
print(" ")

################################################################
#
# RUN PHYLOPHLAN 
#
#################################################################


print(" ")
print(" ")
print('**********************************************')
print(" Step 6: Phlyogenetic inference using \033[92m Phylophlan \033[0m" )
print('**********************************************')
print(" ")
print(" ")


os.system('cp -R /media/data/SANDBOXES/sa01fd/ProkComp/Annotations/' + AnalysisName + ' /opt/softwares/phylophlan/input/Annotations')
os.system('cd /opt/softwares/phylophlan/')
os.system('./phylophlan.py -u Annotations --nproc 16')
os.system('cd /media/data/SANDBOXES/sa01fd/ProkComp/')
os.system('mkdir Phylophlan_output/' + AnalysisName) 
os.system('cp -R /opt/softwares/phylophlan/data/Annotations /media/data/SANDBOXES/sa01fd/ProkComp/Phylophlan_output/' + AnalysisName + '/data')
os.system('cp -R /opt/softwares/phylophlan/output/Annotations /media/data/SANDBOXES/sa01fd/ProkComp/Phylophlan_output/' + AnalysisName + '/output')
os.system('rm -r /opt/softwares/phylophlan/data/Annotations')
os.system('rm -r /opt/softwares/phylophlan/output/Annotations')
os.system('mkdir Phylogenies/' + AnalysisName)
os.system('mkdir ./Phylogenies/'+ AnalysisName + '/alignments')
os.system('cp /media/data/SANDBOXES/sa01fd/ProkComp/Phylophlan_output/' + AnalysisName + '/data/aln.fna ./Phylogenies/' + AnalysisName + '/alignments/phylophlan_alignment.fasta')

################################################################
#
# RAxML on phylophlan 
#
#################################################################

print(" ")
print(" ")
print('**********************************************')
print(" Step 6: ML tree on  phylophlan using \033[92m RAxML \033[0m" )
print('**********************************************')
print(" ")
print(" ")

os.system('mkdir RAxML_output/' + AnalysisName)
os.chdir('./RAxML_output/' + AnalysisName)
os.system('/opt/softwares/RAxML/8.2.9/raxmlHPC-PTHREADS-AVX -f a -s /media/data/SANDBOXES/sa01fd/ProkComp/Phylogenies/' + AnalysisName + '/alignments/phylophlan_alignment.fasta -x 12345 -# 100 -m PROTCATAUTO -n ' + AnalysisName + '_Phylophlan_100_autosubst  -p 12345 -T 16')
os.chdir('/media/data/SANDBOXES/sa01fd/ProkComp/')
os.system('cp ./RAXML_output/' + AnalysisName + '/RAxML_bipartitions.' + AnalysisName + '_Phylophlan_100_autosubst ./Phylogenies/' + AnalysisName + '/trees/phylophlan_RAxML_100')


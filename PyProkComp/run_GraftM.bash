#!/bin/bash

mkdir GraftM_output;
 
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.01.2013_08_greengenes_61_otus.gpkg/ --output_directory ./GraftM_output/test_Marinobacter_16S;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.07.ribosomal_protein_L2_rplB.gpkg/ --output_directory ./GraftM_output/o_4.07_L2_rplB;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.08.ribosomal_protein_L3_rplC.gpkg/ --output_directory ./GraftM_output/o_4.08_L3_rplC;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.09.ribosomal_protein_L5_rplE.gpkg/ --output_directory ./GraftM_output/o_4.09_L5_rplE;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.10.ribosomal_protein_L6_rplF.gpkg/ --output_directory ./GraftM_output/o_4.10_L6_rplF;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.11.ribosomal_protein_L10.gpkg/ --output_directory ./GraftM_output/o_4.11_L10;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.12.ribosomal_protein_L11_rplK.gpkg/ --output_directory ./GraftM_output/o_4.12_L11_rplK;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.13.ribosomal_protein_L14b_L23e_rplN.gpkg/ --output_directory ./GraftM_output/o_4.13_L14b_L23e_rplN;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.14.ribosomal_protein_L16_L10E_rplP.gpkg/ --output_directory ./GraftM_output/o_4.14_L16_L10E_rplP;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.15.ribosomal_protein_S2_rpsB.gpkg/ --output_directory ./GraftM_output/o_4.14_S2_rpsB;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.16.ribosomal_protein_S5.gpkg/ --output_directory ./GraftM_output/o_4.16_S5;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.17.ribosomal_protein_S7.gpkg/ --output_directory ./GraftM_output/o_4.17_S7;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.18.ribosomal_protein_S10_rpsJ.gpkg/ --output_directory ./GraftM_output/o_4.18_S10_rpsJ;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.19.ribosomal_protein_S12_S23.gpkg/ --output_directory ./GraftM_output/o_4.19_S12_S23;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.20.ribosomal_protein_S15P_S13e.gpkg/ --output_directory ./GraftM_output/o_4.20_S15P_S13e;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.21.ribosomal_protein_S19_rpsS.gpkg/ --output_directory ./GraftM_output/o_4.21_S19_rpsS;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.27.methyl_coenzyme_reductase_alpha_subunit.mcrA.gpkg/ --output_directory ./GraftM_output/o_4.27_methyl_coenzyme_reductase_alpha_subunit_mcrA;
graftM graft --forward ./QC_Assemblies/* --graftm_package ./GraftM/data.ace.uq.edu.au/public/graftm/4.40.2013_08_greengenes_97_otus.with_euks.gpkg/ --output_directory ./GraftM_output/o_4.40.2013.8.greengenes;


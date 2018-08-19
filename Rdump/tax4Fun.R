devtools::install_github('EESI/themetagenomics',build_vignettes=TRUE)


source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")



install.packages("themetagenomics", dependencies=TRUE, repos='http://cran.rstudio.com/')


library(themetagenomics)
GEVERS$OTU[1:5,1:5]
GEVERS$TAX[1:5,1:3]

tmp <- tempdir()
download_ref(tmp,reference='gg_ko',overwrite=FALSE)

system.time(FUNCTIONS <- picrust(GEVERS$OTU,rows_are_taxa=FALSE,
                                 reference='gg_ko',reference_path=tmp,
                                 cn_normalize=TRUE,sample_normalize=FALSE,
                                 drop=TRUE))
FUNCTIONS$fxn_table

library(gplots)






testPRISCUIT = read.table('~/DATA/T6SS/TestPRISCUIT.txt',header=TRUE,sep='\t') 
rownames(testPRISCUIT) = testPRISCUIT$X
testPRISCUIT =  testPRISCUIT[,2:length(colnames(testPRISCUIT))]

testPRISCUIT2 = read.table('~/DATA/T6SS/TestPRISCUIT2.txt',header=TRUE,sep='\t') 
rownames(testPRISCUIT2) = testPRISCUIT2$X
testPRISCUIT2 =  testPRISCUIT2[,2:length(colnames(testPRISCUIT2))]


testPRISCUIT
OTUtest = cbind(rep(1,length(rownames(testPRISCUIT))),rep(1,length(rownames(testPRISCUIT))))
rownames(OTUtest) = rownames(testPRISCUIT)
colnames(OTUtest) = c('S1','S2')

system.time(FUNCTIONS <- picrust(BOEVER$OTU,rows_are_taxa=FALSE,
                                 reference='gg_ko',reference_path=tmp,
                                 cn_normalize=TRUE,sample_normalize=FALSE,
                                 drop=TRUE))
BOEVER <- list()                             
BOEVER$TAX = testPRISCUIT
BOEVER$TAX2 = testPRISCUIT2
BOEVER$OTU =  t(OTUtest)                              


download_ref(tmp,reference='silva_ko',overwrite=FALSE)



system.time(FUNCTIONS <- t4f(BOEVER$OTU,rows_are_taxa=FALSE,tax_table=BOEVER$TAX2,
                             reference_path=tmp,type='uproc',short=TRUE,
                             cn_normalize=TRUE,sample_normalize=TRUE,drop=TRUE))

FUNCTIONS$fxn_table
names(FUNCTIONS$fxn_meta)


BOEVER$OTU2 = diag(x = 1, length(rownames(testPRISCUIT2)), length(rownames(testPRISCUIT2)))
rownames(BOEVER$OTU2)=rownames(testPRISCUIT2)
colnames(BOEVER$OTU2)= rownames(testPRISCUIT2)


system.time(FUNCTIONS <- t4f(BOEVER$OTU2,rows_are_taxa=FALSE,tax_table=BOEVER$TAX2,
                             reference_path=tmp,type='uproc',short=TRUE,
                             cn_normalize=TRUE,sample_normalize=TRUE,drop=TRUE))


system.time(FUNCTIONS <- t4f(BOEVER$OTU2,rows_are_taxa=FALSE,tax_table=BOEVER$TAX2,
                             reference_path=tmp,type='uproc',short=TRUE,
                             cn_normalize=TRUE,sample_normalize=TRUE,drop=FALSE))

system.time(FUNCTIONS <- t4f(BOEVER$OTU2,rows_are_taxa=FALSE,tax_table=BOEVER$TAX2,
                             reference_path=tmp,type='pauda',short=FALSE,
                             cn_normalize=FALSE,sample_normalize=FALSE,drop=FALSE))


heatmap.2(FUNCTIONS$fxn_table,scale='none',trace='none',col=colorRampPalette(c("forestgreen", "yellow", "red"))(n = 299))


##### GENES

selectedKO = names(FUNCTIONS$fxn_meta$KEGG_Description[grep('nucleotide',FUNCTIONS$fxn_meta$KEGG_Description)])
selFUNC = FUNCTIONS$fxn_table[, grepl(paste(selectedKO,collapse='|'),colnames(FUNCTIONS$fxn_table))]

heatmap.2(selFUNC , scale='none',trace='none',col=colorRampPalette(c("forestgreen", "yellow", "red"))(n = 299))

##### PATHWAYS

selectedKO = names(FUNCTIONS$fxn_meta$KEGG_Pathways[grep('carbon metabolism',FUNCTIONS$fxn_meta$KEGG_Pathways)])
selFUNC = FUNCTIONS$fxn_table[, grepl(paste(selectedKO,collapse='|'),colnames(FUNCTIONS$fxn_table))]

heatmap.2(selFUNC , scale='none',trace='none',col=colorRampPalette(c('white',"black", "yellow", "red"))(n = 299))


unique(FUNCTIONS$fxn_meta$KEGG_Pathways)



selectedKO = names(FUNCTIONS$fxn_meta$KEGG_Pathways[grep("polysac",FUNCTIONS$fxn_meta$KEGG_Pathways)])
selFUNC = FUNCTIONS$fxn_table[, grepl(paste(selectedKO,collapse='|'),colnames(FUNCTIONS$fxn_table))]

heatmap.2(selFUNC , scale='none',trace='none',col=colorRampPalette(c('white',"black", "yellow", "red"))(n = 299))


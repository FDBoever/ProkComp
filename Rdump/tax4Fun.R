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


heatmap.2(FUNCTIONS$fxn_table,scale='none',trace='none',col=colorRampPalette(c("white","forestgreen", "yellow", "red"))(n = 299))


##### GENES

selectedKO = names(FUNCTIONS$fxn_meta$KEGG_Description[grep('nucleotide',FUNCTIONS$fxn_meta$KEGG_Description)])
selFUNC = FUNCTIONS$fxn_table[, grepl(paste(selectedKO,collapse='|'),colnames(FUNCTIONS$fxn_table))]
rownames(selFUNC)= BOEVER$TAX2[rownames(selFUNC),"Genus"]
heatmap.2(selFUNC , scale='none',trace='none',col=colorRampPalette(c("forestgreen", "yellow", "red"))(n = 299))

##### PATHWAYS

selectedKO = names(FUNCTIONS$fxn_meta$KEGG_Pathways[grep('carbon metabolism',FUNCTIONS$fxn_meta$KEGG_Pathways)])
selFUNC = FUNCTIONS$fxn_table[, grepl(paste(selectedKO,collapse='|'),colnames(FUNCTIONS$fxn_table))]
rownames(selFUNC)= BOEVER$TAX2[rownames(selFUNC),"Genus"]

heatmap.2(selFUNC , scale='none',trace='none',col=colorRampPalette(c('white',"black", "yellow", "red"))(n = 299))


unique(FUNCTIONS$fxn_meta$KEGG_Pathways)



selectedKO = names(FUNCTIONS$fxn_meta$KEGG_Pathways[grep("polysac",FUNCTIONS$fxn_meta$KEGG_Pathways)])
selFUNC = FUNCTIONS$fxn_table[, grepl(paste(selectedKO,collapse='|'),colnames(FUNCTIONS$fxn_table))]
rownames(selFUNC)= as.character(BOEVER$TAX2[rownames(selFUNC),"Genus"])

heatmap.2(selFUNC , scale='none',trace='none',col=colorRampPalette(c('white',"black", "yellow", "red"))(n = 299))


selectedKO = c("K02660","K02659","K06596","K06598","K02658","K02657","K01768","K03651","K10914","K10941","K20968","K13060","K18304","K13061","K18099","K18100","K18101","K12990","K20258","K20259","K01657","K01658","K18000","K18001","K18002","K18003","K20257","K17940","K19735","K07678","K20971","K20972","K07689","K20969","K20970","K03563","K20973","K20974","K20975","K20976","K20977","K20978","K02398","K02405","K11912","K11915","K11890","K11891","K11893","K11913","K11902","K11901","K11900","K11903","K11895","K11907","K13487","K13490","K13488","K13489","K13491","K11444","K21019","K21020","K21021","K21022","K21023","K21024","K20997","K16011","K12992","K20987","K20998","K20999","K21000","K21001","K21002","K21003","K21004","K21005","K21006","K21007","K21008","K21009","K21010","K21011","K21012","K19291","K21025")

selectedKO = c("K20987","K20998","K20999","K21000","K21001","K21002","K21003","K21004","K21005","K21006","K21007","K21008","K21009","K21010","K21011","K21012")

selectedKO = c("K01991","K03328","K05399","K05789","K05790","K06861","K07091","K07257","K07265","K07266","K07271","K08280","K08992","K09688","K09689","K09690","K09691","K09774","K10107","K11719","K11720","K16552","K16554","K16565","K16566","K16567","K16568","K16695","K16696","K16712","K16713","K19363","K19421","K19804","K20920","K20921","K20922","K20946","K20947","K20948","K20949","K20950","K20987","K20988","K20997","K20998","K20999","K21000","K21001","K21002","K21003","K21004","K21005","K21006","K21007","K21008","K21009","K21010","K21011","K21012","K21154","K21556","K21557","K22914","K22915","K22916","K22917","K22921","K22922","K22923","K22924")
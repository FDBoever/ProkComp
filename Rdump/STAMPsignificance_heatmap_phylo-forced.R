
#LOAD MANY STAMP FILES

STAMP1 = read.table(file = '~/adhaerens_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
data.frame(COG=as.character(STAMP1$COG),p=STAMP1$p.values)
STAMP2 = read.table(file = '~/algicola_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP3 = read.table(file = '~/antarcticus_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP4 = read.table(file = '~/excellens_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP5 = read.table(file = '~/hydrocarb_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP6 = read.table(file = '~/psychro_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP7 = read.table(file = '~/vinifirmus_COG_Bonferroni_0.01.tsv', sep = '\t', header = TRUE)
STAMP8 = read.table(file = '~/lipolyticus_enriched_list2.tsv', sep = '\t', header = TRUE)



#load PHYLOGENETIC TREE
RAxMLANVIO = read.tree("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/ANVIO_RAxML/RAxML_bipartitions.autosubst100")
#Root
RAxMLANVIORooted =  midpoint.root(RAxMLANVIO)
#increasing branch lengths
tree2 <- ladderize(RAxMLANVIORooted, right = FALSE)

#can’t have any branch lengths of zero or downstream commands will collapse those nodes…
tree2$edge.length[which(tree2$edge.length == 0)] <- 0.00001
rep_tree_um <- chronopl(tree2,lambda = 0.1,tol = 0)
rep_tree_d <- as.dendrogram(as.hclust.phylo(rep_tree_um))


#	extract accessory genes
AccesoryDF = AnvioData[,!colSums(AnvioData)==1 ]
AccesoryDF  = AccesoryDF[,!colSums(AccesoryDF)==60 ]
RareAccesoryDF  = AccesoryDF[,!colSums(AccesoryDF)<7]

#force row order so that it matches the order of leafs in rep_tree_d
clade_order <- order.dendrogram(rep_tree_d)
clade_name <- labels(rep_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(AccesoryDF))
combined_ordered_matrix <- AccesoryDF[new_order,]


ACCESLAnv = AnvioData[,colSums(AnvioData)>1 & colSums(AnvioData)<60 ]
ACCESmetadata = ANVIO[ANVIO$protein_cluster_id %in% colnames(ACCESLAnv),]
ACCESmetadata2 = ACCESmetadata
ACCESmetadata2 $combined = as.factor(paste(ACCESmetadata2 $COG_FUNCTION, ACCESmetadata2 $COG_FUNCTION_ACC,sep='___'))




######################
# FILTER STAMP DATA
######################
#column 7 are the p-values, cutoff here is 0.01

mat = t(table(droplevels(subset(ACCESmetadata2, COG_FUNCTION_ACC  %in% c(as.character(STAMP1[STAMP1[,7]<0.01,]$COG),as.character(STAMP2[STAMP2[,7]<0.01,]$COG),as.character(STAMP3[STAMP3[,7]<0.01,]$COG),as.character(STAMP4[STAMP4[,7]<0.01,]$COG),as.character(STAMP5[STAMP5[,7]<0.01,]$COG),as.character(STAMP6[STAMP6[,7]<0.01,]$COG),as.character(STAMP7[STAMP7[,7]<0.01,]$COG),as.character(STAMP8[STAMP8[,7]<0.01,]$COG)
))[,c('COG_FUNCTION_ACC','genome_name')])))  

#order the matrix according to the phylogeny
mat_HEAT = mat[,2:ncol(mat)]
mat_HEAT2 = mat_HEAT
#mat_HEAT2[mat_HEAT2>1] <-1
mat_HEAT2_ordered <- mat_HEAT[new_order,]
mat_HEAT2_ordered = as.data.frame.matrix(mat_HEAT2_ordered)


#make dataframe with the significance levels of all groups
annot_df <- data.frame(adhaerens = subset(STAMP1, COG  %in% colnames(mat_HEAT2_ordered))[,7],algicola = subset(STAMP2, COG  %in% colnames(mat_HEAT2_ordered))[,7],antarcticus = subset(STAMP3, COG  %in% colnames(mat_HEAT2_ordered))[,7],excellens = subset(STAMP4, COG  %in% colnames(mat_HEAT2_ordered))[,7],lipo = subset(STAMP8, .  %in% colnames(mat_HEAT2_ordered))[,7],hydrocarb = subset(STAMP5, COG  %in% colnames(mat_HEAT2_ordered))[,7], psychro = subset(STAMP6, COG  %in% colnames(mat_HEAT2_ordered))[,7],vinifirmus = subset(STAMP7, COG  %in% colnames(mat_HEAT2_ordered))[,7])

#fill the dataframe with true/false if they are significant, and change into colors 
annot.df = annot_df>0.01
annot.df[annot.df == TRUE] <- '#808080' 
annot.df[annot.df == FALSE] <- '#FFD700'

#HEATMAP WITH THE SIGNIFICANCE LEVELS ON TOP
library('heatmap.plus')
heatmap.plus(as.matrix(mat_HEAT2_ordered),trace="none",ColSideColors= annot.df
,col = inferno(75),Rowv= rep_tree_d,srtCol=45,cexRow=0.3,cexCol=0.3,main='highest contiruting')

#HEATMAP WITH THE SIGNIFICANCE LEVELS ON TOP
heatmap.2(as.matrix(mat_HEAT2_ordered),trace="none"
,col = inferno(75),Rowv= rep_tree_d,srtCol=45,cexRow=0.3,cexCol=0.3,main='highest contiruting')



###############
# THE ESSENCE
###############

MATRIX = t(table(droplevels(ACCESmetadata2[,c('COG_FUNCTION_ACC','genome_name')])))
reducedMATRIX = MATRIX[,as.character(STAMP2[STAMP2[,7]<0.01,]$COG)]


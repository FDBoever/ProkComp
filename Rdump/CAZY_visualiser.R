CAZY = read.table("~/DATA/MarinobacterGenomics/2018_ProkComp/CAZY_output.txt",header=TRUE)
CAZY2gene = read.table("~/DATA/MarinobacterGenomics/2018_ProkComp/genes.txt",header=TRUE)

CAZY = read.table("~/DATA/MarinobacterGenomics/2018_ProkComp/CAZY_output.txt",header=TRUE)


rownames(CAZY2gene) = CAZY2gene[,1]
CAZY2gene = CAZY2gene[,2:ncol(CAZY2gene)]

rownames(CAZY) = CAZY[,1]
CAZY = CAZY[,2:ncol(CAZY)]

heatmap.2(as.matrix(CAZY),trace='none',col=colorRampPalette(c('white',"#000033", "#FF3300", "#FF3300", "#FF3300"))(n = 50),margins=c(12,12),ColSideColors=ANVIO_cat2$qualityColor)

colnames(CAZY)= gsub("\\.", "_", colnames(CAZY))

CAZYstringent = CAZY[,colnames(CAZY) %in% CheckMANVIO2$Bin_Id]


heatmap.2(as.matrix(CAZYstringent),trace='none',col=colorRampPalette(c('white',"#000033", "#FF3300", "#FF3300", "#FF3300"))(n = 50),margins=c(12,12),ColSideColors=ANVIO_cat2$qualityColor)



clade_order <- order.dendrogram(pruned_tree_d)
clade_name <- labels(pruned_tree_d)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, colnames(CAZYstringent))
combined_ordered_matrix <- CAZYstringent[, new_order]

heatmap.2(as.matrix((combined_ordered_matrix)),trace="none"
,col=colorRampPalette(c('white',"#000033", "#FF3300", "#FF3300", "#FF3300"))(n = 50),srtCol=45,cexRow=0.4,cexCol=0.6,main='',margins=c(10,20),ColSideColors= as.character(cladecolorsHeat[colnames(combined_ordered_matrix)]),Colv= pruned_tree_d)

combined_ordered_matrix[rowSums(combined_ordered_matrix)>5,]
heatmap.2(as.matrix((combined_ordered_matrix[rowSums(combined_ordered_matrix)>15,])),trace="none"
,col=colorRampPalette(c('white',"#000033", "#FF3300", "#FF3300", "#FF3300"))(n = 50),srtCol=45,cexRow=0.4,cexCol=0.6,main='',margins=c(10,20),ColSideColors= as.character(cladecolorsHeat[colnames(combined_ordered_matrix)]),Colv= pruned_tree_d)

CAZY2gene[rownames(combined_ordered_matrix),]


combined_ordered_matrix2 = combined_ordered_matrix
rownames(combined_ordered_matrix2)=paste(rownames(combined_ordered_matrix2),CAZY2gene[rownames(combined_ordered_matrix),'hmmsource'],sep= ' - ')


heatmap.2(as.matrix((combined_ordered_matrix2)),trace="none"
,col=colorRampPalette(c('white',"#000033", "#FF3300", "#FF3300", "#FF3300"))(n = 50),srtCol=45,cexRow=0.4,cexCol=0.6,main='',margins=c(10,20),ColSideColors= as.character(cladecolorsHeat[colnames(combined_ordered_matrix)]),Colv= pruned_tree_d)


cazyGroups = c()
for( i in c('GH','GT','PL','CE','CBM')){
	cazyGroups = cbind(cazyGroups , colSums(combined_ordered_matrix[grepl(i,rownames(combined_ordered_matrix)),]))
}
colnames(cazyGroups) = c('GH','GT','PL','CE','CBM')

heatmap.2(as.matrix(t(cazyGroups)),trace="none", col=colorRampPalette(c('white',"#000033", "#FF3300"))(n = 50),srtCol=45,cexRow=1.5,cexCol=0.6,main='',margins=c(10,20),ColSideColors= as.character(cladecolorsHeat[colnames(combined_ordered_matrix)]),Colv= pruned_tree_d)

heatmap.2(as.matrix(t(cazyGroups)),trace="none", scale=c('row')
,col=colorRampPalette(c("#000033", 'white',"#FF3300"))(n = 50),srtCol=45,cexRow=1.5,cexCol=0.6,main='',margins=c(10,20),ColSideColors= as.character(cladecolorsHeat[colnames(combined_ordered_matrix)]),Colv= pruned_tree_d)

cazyCAT = ANVIO_cat2
rownames(cazyCAT) = cazyCAT[,1]
cazyCAT = cazyCAT[,2:ncol(cazyCAT)]



cazyGroupsPhylo = cbind(cazyGroups ,cazyCAT[rownames(cazyGroups),])


moltenCAZY = melt(cazyGroupsPhylo[,c('GH','GT','PL','CE','CBM','group2')],id.var = ('group2'))


ggboxplot(moltenCAZY, x = "group2",color="group2",
          y = c("value"),
          combine = TRUE,
          add = "jitter",                              # Add jittered points
          add.params = list(size = 0.1, jitter = 0.2)  # Point size and the amount of jittering
          ) + scale_colour_manual(values = c("#39811D","#86DA81","#A6D8D4","#F2A968","#F2EC70","#E38FDD","#898989","#76AECF","#B34D22"))  +stat_compare_means() + facet_wrap(~variable,scales='free')
          
          
my_comparisons <- list( c("algicola", "antarcticus"), c("algicola", "hydrocarbo"), c("algicola", "psychro") )

p<-ggplot(moltenCAZY,aes(group2, value,colour= group2))+ylab("nr of CAZYmes")
p<-p+geom_boxplot()+geom_jitter()+theme_bw()+
 facet_wrap( ~ variable , scales='free',nrow=1)
p<-p+ scale_colour_manual(values = c("#39811D","#86DA81","#A6D8D4","#F2A968","#F2EC70","#E38FDD","#898989","#76AECF","#B34D22"))+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90))+stat_compare_means(comparisons = my_comparisons)+stat_compare_means()
print(p)






library(ggplot2)
df = read.table('~/DATA/T6SS/MIC.txt',header=TRUE,sep='\t')
df$diff = df$t2-df$t1



for(ab in unique(df$AB)){
	df2 = df[df$AB==ab,]
	g = ggplot(df2, aes(y = factor(Row, rev(levels(Row))),x = factor(Col))) + 
     geom_point(aes(fill = diff),colour="black", size =18, shape=21)  +theme_bw() +
     labs(x=NULL, y = NULL)+scale_fill_gradient2(low = "firebrick", mid = "white", high = "forestgreen", midpoint = 0)+ggtitle(ab)+theme(axis.text=element_text(size=21),
        axis.title=element_text(size=24))
	pdf(paste("~/DATA/T6SS/plot_MIC_", ab,".pdf",sep=""),width=10, height=6)
  	print(g)
  	dev.off()

}




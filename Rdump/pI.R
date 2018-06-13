library(ggplot2)
library(plyr)



#Read file
pI = read.table("~/Genomics/ProkComp/pI.txt",header=TRUE)
pI = read.table("~/Genomics/ProkComp/pIbacterioplankton.txt",header=TRUE)

##################################################
######## VALUABLE FUNCTIONS #########
# multiplot() will allow you to visualize multiple plots in a grid arrangment

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##################################################



#add additional column to classify protein as basic vs acidic
#play with the threshold number
pI$acidic <- ifelse(pI$pI > 7.5 ,"basic", "acidic")


#for up to 3 proteomes
#g2 = ggplot(pI, aes(x= pI, color= genome, fill= genome)) +
#geom_histogram(aes(y=..density..), position="identity", alpha=0.2)+
#geom_density(alpha=0.2)+
#scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
#scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
#labs(title="",x="pI", y = "Density")+
#theme_classic()

#global view, shows all the proteomes on one graph
g1 = ggplot(pI, aes(x= pI, color= "black", fill= "black")) +
geom_histogram(aes(y=..density..), position="identity", alpha=1,binwidth = 0.1)+ scale_y_continuous(expand = c(0, 0))+scale_color_manual(values=c("black")) + scale_fill_manual(values=c("black"))+ 
labs(title= strBias ,x="class of pI", y = "Density")+
theme_classic()+theme(legend.position="none")


g2 = ggplot(pI, aes(x=pI,y= log(length),color = acidic)) + geom_point()+
scale_color_manual(values=c("grey", "black")) + scale_fill_manual(values=c("grey", "black"))+ labs(title=" ",x="pI", y = "log Length")+ theme_classic()+theme(legend.position="none")


multiplot(g1,g2,cols=2)



##################################################
# plot visialisation per proteome


myplots <- list()
myplots2 <- list()
geno= c()
bias = c()

i = 1
for (organism in unique(pI$genome)){
	pIsub =subset(pI, genome ==organism)
	Nbasic = length(which(pIsub$acidic == 'basic'))
	Nacidic = length(which(pIsub$acidic == 'acidic'))
	pIBias = (Nbasic-Nacidic)/(Nbasic+Nacidic)*100
	strBias = paste(organism,":",round(pIBias,2),"%",sep=" ")
	
	bias = c(bias, pIBias)
	geno = c(geno, organism)
	p1 <- ggplot(pIsub, aes(x= pI, color= "black", fill= "black")) + geom_histogram(aes(y=..density..), position="identity", alpha=1,binwidth = 0.1) + scale_y_continuous(limits = c(0,0.6), expand = c(0, 0))+scale_color_manual(values=c("black")) + scale_fill_manual(values=c("black"))+ labs(title= strBias, x="class of pI", y = "Density")+ theme_classic()+theme(legend.position="none",plot.title = element_text(size=6))
	g2 = ggplot(pIsub, aes(x=pI,y= log(length),color = acidic)) + geom_point()+scale_color_manual(values=c("grey", "black")) + scale_fill_manual(values=c("grey", "black"))+ labs(title= strBias,x="pI", y = "log Length")+ theme_classic()+theme(legend.position="none",plot.title = element_text(size=6))
	myplots[[i]] <- p1 
	myplots2[[i]] <- g2 
	i = i + 1
}
#stores all the calculated b-values in a new data.frame
bias_data = data.frame(cbind(genome=geno,b=bias))


#play arround howmany of the plots you want 
multiplot(plotlist = myplots[1:12],cols=6)
multiplot(plotlist = myplots,cols=9)
multiplot(plotlist = myplots2,cols=9)
multiplot(plotlist = myplots,cols=6)
multiplot(plotlist = myplots,cols=9)


##################################################
# PI BIASS per genome

q = ggplot(bias_data, aes(x=reorder(genome,bias), y=bias)) + 
    geom_point() +      # Thinner lines
    xlab("Genome") +
    ylab("pI bias") +
    ggtitle("pI bias") +
    theme_bw()
q + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5)) + coord_cartesian(ylim = c(45, 70)) + coord_flip()

#show violin graph
 p <- ggplot(pI, aes(factor(genome), pI))
p + geom_violin()+theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5)) + coord_cartesian(ylim = c(45, 70)) + coord_flip()


##################################################

# Select genomes of interest from the set
pIsub1 =subset(pI, genome == "Marinobacter_lipolyticus_SM19.faa")
pIsub2 =subset(pI, genome == "Marinobacter_lipoliticus_BF04_CF-4.faa")
rbind(pIsub1,pIsub2)


g2 = ggplot(rbind(pIsub1,pIsub2), aes(x= pI, color= genome, fill= genome)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.2,binwidth=0.1)+
geom_density(alpha=0.2)+scale_y_continuous(limits = c(0,0.6), expand = c(0, 0))+
scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
labs(title="",x="pI", y = "Density")+
theme_classic()

g3 = ggplot(rbind(pIsub1,pIsub2), aes(x= pI, color= genome, fill= genome)) +
geom_density(alpha=0.2)+scale_y_continuous(limits = c(0,0.6), expand = c(0, 0))+
scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
labs(title="",x="pI", y = "Density")+
theme_classic()

multiplot(g2,g3,cols=2)


pI = read.table("~/Genomics/ProkComp/AA_composition.txt",header=TRUE)





data1=data.frame(cbind(pI[1],pI[7:25]))
head(data1)

longdf=melt(data1,id.vars="genome")
head(longdf)

df=ddply(longdf,.(genome,variable),summarize,mean(value))
colnames(df)[length(df)]="value"
df

ggplot(data=df,aes(x=variable,y=value,group=genome,colour=genome))+
        geom_point()+geom_line()

ggplot(data=df,aes(x=variable,y=value,group=genome,colour=genome,fill=genome))+
        geom_point()+geom_polygon(alpha=0.001)+coord_polar()

ggplot(data=df,aes(x=variable,y=value,group=genome))+
        geom_point()+geom_polygon(alpha=0.001)+coord_polar()



ggplot(data=df,aes(x=variable,y=value,group=genome,colour=genome))+
        geom_point()+geom_polygon(alpha=0.1)+coord_polar()


rescale_df=function(data,groupvar=NULL){
        if(is.null(groupvar)) df=data
        else df=data[,-which(names(data) %in% groupvar)]
        
        select=sapply(df,is.numeric)
        df[select]=lapply(df[select], scales::rescale)
        if(!is.null(groupvar)) {
                df=cbind(df,data[[groupvar]])
                colnames(df)[length(df)]=groupvar
        }        
        df
}

rescaled=rescale_df(data1)
head(rescaled)


longdf2=melt(rescaled,id.vars="genome")
head(longdf2)


df2=ddply(longdf2,.(genome,variable),summarize,mean(value))
colnames(df2)[length(df2)]="value"
df2


ggplot(data=df2,aes(x=variable,y=value,group=genome,colour=genome,fill=genome))+
        geom_point()+geom_polygon(alpha=0.4)+coord_polar()


coord_radar <- function (theta = "x", start = 0, direction = 1) 
{
        theta <- match.arg(theta, c("x", "y"))
        r <- if (theta == "x") 
                "y"
        else "x"
        ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
                direction = sign(direction),
                is_linear = function(coord) TRUE)
}



ggplot(data=df2,aes(x=variable,y=value,group=genome,colour=genome,fill=genome))+
        geom_point()+geom_polygon(alpha=0.4)+coord_radar()


ggplot(data=df2,aes(x=variable,y=value,group=genome,colour=genome,fill=genome))+
        geom_point()+geom_polygon(alpha=0.4)+coord_radar()+ylim(0,1)+
        theme(legend.position="bottom")+xlab("")+ylab("")

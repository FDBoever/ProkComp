install.packages('xlsReadWrite',dependencies=TRUE)
install.packages('plyr',dependencies=TRUE)
library(reshape)

require(xlsx)
read.xlsx("myfile.xlsx", sheetName = "Sheet1")
library(gdata)
library(xlsx)
library(plyr)


path = "C:/Users/SA03NB/Documents/Nephelo/Tesr_R/Gradientest_plasticplate_n_6/55_10_platen_6/"
out.file<-""
file.names <- dir(path, pattern =".xlsx")

DF = cbind()

for(i in 2:length(file.names)){
 	file = read_excel(paste(path,file.names[i],sep=''), skip = 11,col_names=TRUE)
	file = data.frame(file)	
	print(file)
	DF = cbind(DF, file[,4])
}

colnames(DF) = file.names[2:length(file.names)]
DF = cbind( file[,1:3], DF)

dilution  = c(rep(c(75,50,25,12.5,6,3),3),rep('blank',3),rep(2,3))

DF = cbind(dilution  , DF)





Blanks = DF[DF$dilution=='blank',5:dim(DF)[2]]
subtractedDF = DF[DF$dilution!='blank',5:dim(DF)[2]]
meanBlanks = colSums(Blanks)/3

#subtract blanks means per column (file)
subtractedDF  = subtractedDF  - meanBlanks

#subtract blanks overal averaged for all files
#subtractedDF  = subtractedDF - mean(meanBlanks)

dilution2  = c(rep(c(75,50,25,12.5,6,3),3),rep(2,3))
subtractedDF  = cbind(dilution2 , subtractedDF  )


full2 = melt(subtractedDF ,id=c("dilution2"))

#overall means and sd
ds <- plyr::ddply(full2 , c("dilution2"),plyr::summarise, mean = mean(value), sd = sd(value))

#mean and sd per file, column
ds <- plyr::ddply(full2 , c("dilution2",'variable'),plyr::summarise, mean = mean(value), sd = sd(value))


ggplot(ds,las=3,aes(x=dilution2,y=mean, group= variable)) + ylab('nephelometer') + xlab("dilution") + geom_point(aes(color= variable)) + geom_line(aes(color= variable)) + geom_errorbar(aes(ymin = mean-sd,ymax = mean +sd , color = variable),width = 0.01)  + scale_color_grey() + theme_classic() + theme(legend.position="none")

ggplot(full2,aes(x=dilution2,y=value)) + ylab('nephelometer') + xlab("dilution") + geom_point(aes(color= variable)) + geom_smooth() + theme_classic() + theme(legend.position="none")



test.mod1 = lm(value~ dilution2, data = full2)
summary(test.mod1)


Rsquared = summary(test.mod1)$adj.r.squared



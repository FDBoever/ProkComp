setwd("/Users/FrederikDeBoever")
full = read.table("~/Desktop/Desktop backup/FDB_thesis 3/Data as Textfiles/Screening_Wide.txt",header=FALSE)

install.packages('grofit',dependencies=T)
library(grofit)
#create dataframe suitable for grofit analysis.

unique(full$V11)

timeseries = c(0,12,36,60,84,108,132)

time = matrix(rep(timeseries,180), 180, 7, byrow=T)
data = 0
data = data.frame(
exp_id=paste(full$V1),
info=rep(c("rep1","rep2", "rep3","rep4"), 45),
conc=rep(1,180)  # (I set the concentration value to an arbitrary value of 1)
)

tpoints = t(rbind(full$V3,full$V4,full$V5,full$V6,full$V7,full$V8,full$V9,full$V10))
data = cbind(data, tpoints)

#grofit
library(grofit)


for (j in 1:(nrow(tpoints)/4) ){
	plot (time, tpoints, pch=20, cex=0.8, col="gray", ylim=c(0, 0.3))
	for (i in 1:4 ) {
	gr = gcFitSpline(timeseries,tpoints[i*j,])
	points(gr$fit.time,gr$fit.data,t='l')
	}
}




par(mfrow=c(5,9),mar=c(c(0, 0, 0, 0)),oma = c(4.5, 4.5, 1.5, 1.5),tcl = -0.25)

for (j in 1:(nrow(tpoints)/4) ){
	if(j %in% seq(1, 36, 9)){
	plot (time[1:4,], tpoints[(j*4-3):(4*j),], pch=20, cex=0.8, col="darkorange3", ylim=c(0, 0.5),xaxt="n")
	}
	else if(j %in% c(38:45)){
	plot (time[1:4,], tpoints[(j*4-3):(4*j),], pch=20, cex=0.8, col="darkorange3", ylim=c(0, 0.5),yaxt="n")
	}
	else if(j ==37){
	plot (time[1:4,], tpoints[(j*4-3):(4*j),], pch=20, cex=0.8, col="darkorange3", ylim=c(0, 0.5))
	}
	else{
	plot (time[1:4,], tpoints[(j*4-3):(4*j),], pch=20, cex=0.8, col="darkorange3", ylim=c(0, 0.5),xaxt="n",yaxt="n")
	}
	if(j!=100){
	for (i in 1:4 ) {
	gr = gcFitModel(timeseries,tpoints[i+(4*j)-4,])
	points(gr$fit.time,gr$fit.data,t='l',col='darkorange3')
	axo =  gcFitModel(timeseries, tpoints[i,])
	points(axo$fit.time,axo$fit.data,t="l",col="gray")
	points(axo$fit.time,axo$fit.data,col="gray",pch=20, cex=0.8)
	#gr2 = gcFitSpline(timeseries,tpoints[i+(4*j)-4,])
	#points(gr2$fit.time,gr2$fit.data,t='l',col="gray38")
	}}
	else {for (i in 1:4 ) {
		gr2 = gcFitSpline(timeseries,tpoints[i+(4*j)-4,])
	points(gr2$fit.time,gr2$fit.data,t='l',col="gray38")
	}}
	text(40,0.46,data[j*4,1])

}
mtext('Time (h)', side = 1, outer = TRUE, line = 2)
mtext('minimal fluorescence (f0)', side = 2, outer = TRUE, line = 2)



#EXTRACT THE DATA

non_interactive = grofit.control(interactive=FALSE)
TestRun = grofit(time,data, FALSE, non_interactive)

# BASED ON SPLINE
mu.spline <- TestRun$gcFit$gcTable$mu.spline
lambda.spline <- TestRun $gcFit$gcTable$lambda.spline
Assymp.spline <- TestRun $gcFit$gcTable$A.spline
mu.model <- TestRun$gcFit$gcTable$mu.model
lambda.model <- TestRun $gcFit$gcTable$lambda.model
Assymp.model <- TestRun $gcFit$gcTable$A.model
full = cbind(full,mu.spline,Assymp.spline,lambda.spline,mu.model,Assymp.model,lambda.model)



#get other data

timeshort = time[1,1:4]
datashort = tpoints[1,1:4]
testmatrix =data.frame(timeshort,datashort)

Mod.B =c()
modList = c()
maxf0 = c()
meanf0 =c()

for (j in 1:(nrow(tpoints)) ){
	timeshort = time[j,1:4]
	datashort = tpoints[j,1:4]
		datashort2 = tpoints[j,5:7]

	testmatrix =data.frame(timeshort,datashort)
	mod <- nls(datashort ~ exp(a + b * timeshort), data = testmatrix, start = list(a = 0, b = 0))
	Mod.B =  c(Mod.B,coef(mod)["b"])
	#modList = c(modList, Mod.B)
	maxf0 = c(maxf0,max(datashort2))
	meanf0 = c(meanf0,mean(datashort2))
}


full = cbind(full,Mod.B)
full = cbind(full,maxf0,meanf0)
#plot graphs 




library(lattice)
par(mfrow = c(2, 1), mar = c(4.1, 4.1, 2.1, 1.1))

plot(full$mu.model ~ full$V1,las=2,cex.axis =0.8,cex.lab=0.8,xlab='treatment',ylab='mu')
plot(full$Assymp.model ~ full$V1,las=2,cex.axis =0.8,cex.lab=0.8,xlab='treatment',ylab='Assymptote')
#bwplot(data$lambda ~ data$exp_id,las=2)

par(mfrow = c(2, 1), mar = c(4.1, 4.1, 2.1, 1.1))

plot(full$Mod.B ~ full$V1,las=2,cex.axis =0.8,cex.lab=0.8,xlab='treatment',ylab='mu')
plot(full$maxf0 ~ full$V1,las=2,cex.axis =0.8,cex.lab=0.8,xlab='treatment',ylab='Assymptote')
#bwplot(data$lambda ~ data$exp_id,las=2)

______________
#ORDERED




full3 =subset(full, V1!="yeast")
full3$V1 <- factor(full3$V1)

library(lattice)
par(mfrow = c(2, 1), mar = c(4.1, 4.1, 2.1, 1.1))
plot(full3$mu.model ~ reorder(full3$V1,full3$mu.model),las=2,cex.axis =0.8,cex.lab=0.8,xlab=' ',ylab='growth rate (grofit)')
plot(full3$Assymp.model ~ reorder(full3$V1,full3$Assymp.model),las=2,cex.axis =0.8,cex.lab=0.8,xlab=' ',ylab='Assymptote (grofit)')
#bwplot(data$lambda ~ data$exp_id,las=2)

par(mfrow = c(2, 1), mar = c(4.1, 4.1, 2.1, 1.1))

plot(full3$Mod.B ~ reorder(full3$V1,full3$Mod.B),las=2,cex.axis =0.8,cex.lab=0.8,xlab='',ylab='growth rate (nmle)')
plot(full3$maxf0 ~ reorder(full3$V1,full3$maxf0),las=2,cex.axis =0.8,cex.lab=0.8,xlab='',ylab='Assymptote (maxf0)')
#bwplot(data$lambda ~ data$exp_id,las=2)



________________




plot(sort(full$mu.model) ~ full$V1,las=2,cex.axis =0.8,cex.lab=0.8,xlab='treatment',ylab='mu')


par(mfrow = c(2, 2), mar = c(4.1, 4.1, 2.1, 1.1))
plot(full$mu.model,full$mu.spline,cex.axis =0.7,cex=0.8,pch=20,xlab="model(mu)",ylab="spline(mu)",cex.lab=0.9)
mod1 = lm(formula =  full$mu.spline ~ full$mu.model)
abline(mod1)
text(0.0027,0.008,paste("Rsq=",round(summary(mod1)$r.squared,4)),cex=0.7)
plot(full$Assymp.model,full$Assymp.spline,cex.axis =0.7,cex=0.8,pch=20,,xlab="model(A)",ylab="spline(A)",cex.lab=0.9)
mod1 = lm(formula =  full$Assymp.spline ~ full$Assymp.model)
abline(mod1)
text(0.18,0.42,paste("Rsq=",round(summary(mod1)$r.squared,4)),cex=0.7)





#par(mfrow = c(1, 2), mar = c(4.1, 4.1, 2.1, 1.1))

plot(full$mu.model,full$Mod.B,xlab="mu:nlme()",ylab="mu:grofit(model)",cex.axis =0.7,cex=0.8,pch=20,cex.lab=0.9)
mod1 = lm(formula =  full$Mod.B ~ full$mu.model)
abline(mod1)
text(0.003,0.029,paste("Rsq=",round(summary(mod1)$r.squared,4)),cex=0.7)

plot(full$mu.spline,full$Mod.B,xlab="mu:nlme()",ylab="mu:grofit(spline)",cex.axis =0.7,cex=0.8,pch=20,cex.lab=0.9)
mod1 = lm(formula =  full$Mod.B ~full$mu.spline)
abline(mod1)
text(0.0027,0.029,paste("Rsq=",round(summary(mod1)$r.squared,4)),cex=0.7)








#get means and SD
aggregate(full,by=list(full$V1),function(x) c(mean = mean(x), sd = sd(x)))


________________


f2 = aggregate(mu.model ~ V1, data = full, mean)
f3 = f2[ order(f2$mu.model,f2$V1), ]

x = 


full2= full[ order(full$mu.model, full$Assymp.model, full$V1), ]


boxplot(full$mu.model~reorder(full$V1,full$mu.model),data=df)




# plot sorted data

tagg = aggregate(maxf0 ~ V1, data = full, mean)
t5g = aggregate(maxf0 ~ V1, data = full, sd)

sd = t5g$maxf0
sortinglist = cbind(tagg,x)

sortedlist = sortinglist[ order(sortinglist $maxf0, sortinglist $V1, sortinglist$x), ]

sortlist = read.table("sortedList_exp1.txt",header=FALSE)
plot(sortlist$V2 ~ factor(sortlist$V1),las=2,cex.axis =0.8,cex.lab=0.8,xlab='treatment',ylab='mu')

sortlist

#tagg2 = tagg[ order(tagg$maxf0, tagg$V1), ]
#t5gg2 = t5g[ order(t5g$maxf0, t5g$V1), ]

#plot(tagg2$maxf0 ~ tagg2$V1,las=2,cex.axis =0.8,cex.lab=0.8,xlab='treatment',ylab='mu')



______________________




anov = aov(full$Mod.B ~ full$V1)
TukeyHSD(anov)
plot(full$Mod.B  ~ full$V1,las=2,cex.axis =0.8,cex.lab=0.8,xlab=' ',ylab='b (nmle)')
summary(anov)



anov = aov(full$maxf0 ~ full$V1)
TukeyHSD(anov)
plot(full$maxf0  ~ full$V1,las=2,cex.axis =0.8,cex.lab=0.8,xlab=' ',ylab='max.f0')
summary(anov)

plot(full$lambda  ~ full$V1,las=2,cex.axis =0.8,cex.lab=0.8,xlab=' ',ylab='max.f0')






library(multcomp)
anova = aov(formula = Assymp.spline~ V1, data=full)
summary(glht(anova, linfct=mcp(V1="Dunnett")))


___________________
growth rates on shorter time intervall
___________________


#get other data

timeshorter = time[1,1:3]
datashorter = tpoints[1,1:3]
testmatrixshort =data.frame(timeshorter,datashorter)

Mod.B =c()
modList = c()
maxf0 = c()
meanf0 =c()

for (j in 1:(nrow(tpoints)) ){
	timeshort = time[j,1:2]
	datashort = tpoints[j,1:2]
		datashort2 = tpoints[j,5:7]

	testmatrix =data.frame(timeshort,datashort)
	mod <- nls(datashort ~ exp(a + b * timeshort), data = testmatrix, start = list(a = 0, b = 0))
	Mod.B =  c(Mod.B,coef(mod)["b"])
	#modList = c(modList, Mod.B)
	maxf0 = c(maxf0,max(datashort2))
	meanf0 = c(meanf0,mean(datashort2))
}


full = cbind(full,Mod.B)
full = cbind(full,maxf0,meanf0)
#plot graphs 

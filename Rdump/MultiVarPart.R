# Frederik De Boever
# Last edited: 22/05/18

################################################################################################################################

#	MultiVarPart.R

# This script aims to analysis and explore multiple related datasets, like species abundance and environmental data (and more)  
# 
#	
#
#	CCA 		Canonical correspondence analysis
#			Function cca performs correspondence analysis, or optionally constrained correspondence analysis (a.k.a. canonical correspondence analysis), or optionally partial constrained correspondence analysis
#
#	adonis 	Permutational Multivariate analysis of Variance using distance matrices
#			Analysis of variance using distance matrices — for partitioning distance matrices among sources of variation and fitting linear models (e.g., factors, polynomial regression) to distance matrices; 
#			uses a permutation test with pseudo-F ratios. 
#
#	varpart	variation partitioning
#			The function partitions the variation in community data or community dissimilarities with respect to two, three, or four explanatory tables, using adjusted R−squared in redundancy analysis ordination (RDA) or 			distance-based redundancy analysis. 
#			If response is a single vector, partitioning is by partial regression. Collinear variables in the explanatory tables do NOT have to be removed prior to partitioning.
#

# traits:	 is a generic name for the dataframe containing either; 
#		species abundance, 
#		gene abundance, 
#		phenotypic trait data (such as morphology), 
#		genotypic data (allellic composition) 
#		etc....
# 
# environment:	contains the meta_data information we want to use for CCA
#		Environmental parameters: a classical CCA in ecology

# references:
# 	http://userweb.eng.gla.ac.uk/umer.ijaz
#	http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html
#	https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/varpart
#	https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/cca

# GOOD READS to understand the basics and INTERPRETATION of CCA
#	http://cc.oulu.fi/~jarioksa/opetus/metodi/ordination101.html#65
#	ter Braak, C.J.F. & Verdonschot, P.F.M. Aquatic Science (1995) 57: 255. https://doi.org/10.1007/BF00877430
#	http://ordination.okstate.edu/varpar.html

################################################################################################################################

library(vegan)
library(ggvegan)
library(grid)


traits<-read.csv("~/Downloads/SPE_pitlatrine.csv",row.names=1,check.names=FALSE)

# Transpose the data to have sample names on rows
 traits<-t(traits)

environment<-read.csv("~/Downloads/ENV_pitlatrine.csv",row.names=1,check.names=FALSE)

# ensure that the samples in environment are in the same order as in traits
environment<-environment[rownames(traits),]
# Filter out any samples taxas that have zero entries 

traits<-subset(traits,rowSums(traits)!=0)
# !!!! OPTIONAL !!!!! (need to understand the difference)
# Convert to relative frequencies if working with abundance data
# may want to scale the data here if working with morphometric data?
# traits<-traits/rowSums(traits)


################################################################################################################################
#
# ADONIS : Permutational Multivariate analysis of Variance using distance matrices
#
################################################################################################################################
# as part of the CCA analysis, one could use adonis to select those environmental parameters that are found significant (I personally would not recoomend this)

# Use adonis to find significant environmental variables
traits.adonis <- adonis(traits ~ ., data=environment)

#Extract the "best" variables using a cut-off of interest!
Pr.cutoff = 0.01
bestEnvVariables<-rownames(traits.adonis$aov.tab)[traits.adonis$aov.tab$"Pr(>F)"<=0.01]

#Last two are NA entries, so we have to remove them
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]


################################################################################################################################
#
# CCA: Canonical correspondence analysis
#
################################################################################################################################


# MAKE YOUR CHOISE HERE! Determine the environmental parameters you want to include!

# EITHER WE LOOK AT the complete datasets (including all environmental parameters) (what I'd recommend)
cca_obj<-cca(traits ~ ., data=environment)

# OR 	1) use only those environmental variables in cca that were found significant (BASED ON ADONIS)
# eval(parse(text=paste("cca_obj <- cca(traits ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=environment)",sep="")))

# OR 	2) come up with a subset of environmental parameters which are not co-linear (by using for example VIF (variance inflation factors)-based approaches)
# cca_obj<-cca(traits ~ pH+Temp+TS, data=environment)


#Informative is to compare unconstrained versus constraint (environment) analysis
#--------------------------------------------
#here we use the full environmental dataset as an example
#Compare the output of these two ordinations. In particular, see how the inertia and and rank (number of axes) are decomposed in constrained ordination. 

#unconstrained
m <- cca(traits)

#unconstrained
#mm <- cca(traits, environment)
mm <-cca(traits ~ ., data=environment)

m
mm

#It is often a bad idea to use constrained ordination if you have a very large number of constraints: probably you will not constrain at all, but your analysis may be very similar to ordinary unconstrained analysis that could have been used just as well. You can inspect this with Procrustes rotation: 
plot(procrustes(m, mm))

#You can also see how similar the ordinations are by fitting environmental variables to unconstrained ordination. 
plot(envfit(m, environment))

#One problem with model building is that constraining variables are not in- dependent, but they are correlated. Any correlated variable can be explained with other variables. Such variables are redundant (“expendable”) when they are with other variables, but they may be the best variables alone and prevent other variables from entering the model. A statistic describing this is called variance inflation factor (VIF) which is 1 for completely independent variables, and values above 10 or 20 (depending on your taste) are regarded as highly multicollinear (dependent on others). The VIF of a variable will depend on the set it sits with: 
#NOTE more recent publications on VIF, suggests lower treshholds, (<5), I'd suggest to incorporate some thought process as well, to leave in those that rationally make sense.


VIF.CCA = vif.cca(mm)

#If you want to extract those below a threshold use this
vif.threshold = 5

VIF.CCA[VIF.CCA <vif.threshold]
names(VIF.CCA[VIF.CCA <vif.threshold])


#Package vegan has several alternative types of significance tests. They all can be performed with a function called anova. The name is somewhat mis- leading: the test are based on permutations although the layout of the results is similar as in the standard ANOVA table. The default is an overall test of all variables together: 
anova(mm, permu=200)
anova(mm, by="term", permu=200)


#inspect the perutation tests graphically

df_permutation = attr(anova(mm, by="term", permu=200),'F.perm')[,]
colnames(df_permutation ) =  attr(mm$terms, "term.labels")
molten_permutation = melt(df_permutation)
colnames(molten_permutation) = c('perm','variable','value')

ggplot(molten_permutation, aes(value, fill = variable)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.4) +facet_wrap(~ variable,ncol=1)

ggplot(molten_permutation, aes(value, fill = variable)) + 
  stat_density(aes(y = ..density..), position = "identity", color = "black", alpha = 0.4)


#		!!!!!!!!!! CAUTION !!!!!!!!!!!!!
# make sure you are confident in chosing the environmental parameters to include, typically avoid co-linearity
# once you have your set of variables, rerun the cca and continue

#AFTER YOU MADE YOUR DESCISION WHICH VARIABLES TO KEEP, INCLUDE IN CCA, AND CONTINUE HERE
#--------------------------------------------------------------------------------------------

#standard visualisation as part of vegan package
plot(cca_obj)
plot(cca_obj, display = c("lc","bp"))


#looking at the object to derive statistical metrics
# have a look here for details on how to interpret the inertia https://stackoverflow.com/questions/22535943/how-to-interpret-cca-vegan-output
summary(cca_obj)

#storing the output in an object used for ggplot visualisation
scrs<-scores(cca_obj,display=c("sp","wa","lc","bp","cn"))
attributes(scrs)


# ALTERVATIVE VISUALISATION
#--------------------------------------------

#	ggvegan
#----------------
library(ggvegan)

autoplot(cca_obj)


# custom ggplot visualisation
#---------------------------

# we make a dataframe based on the CCA1 and CCA2 values, but want to annotate these to colorate the dots
# Here make sure you tweak this around, as it is made to work with the testing data I used
# @KATI, one can extract the CCA1 and CCA2 values for plotting, but try and squeeze columns that include DEPTH and TIME?

#Extract site data first
df_sites<-data.frame(scrs$sites,t(as.data.frame(strsplit(rownames(scrs$sites),"_"))))
colnames(df_sites)<-c("x","y","Country","Latrine","Depth")

#Draw sites
p<-ggplot() + geom_point(data=df_sites,aes(x,y,colour=Country))

multiplier <- vegan:::ordiArrowMul(scrs$biplot)

df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)
p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                 arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)
p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

#Either choose text or points
#p<-p+geom_text(data=df_species,aes(x,y,label=rownames(df_species)))
#p<-p+geom_point(data=df_species,aes(x,y,shape="Species"))+scale_shape_manual("",values=2)

p<-p+theme_bw()
print(p)


################################################################################################################################
#
# Variance partitioning 
#
################################################################################################################################
#The basic partial ordination gives decomposition of the inertia into conditional (partialled out) and (residual) constrained component. The decomposition is se- quential: conditions are partialled out before analysing constraints. Reversing the order of terms changes the components. Often we want to have a neutral de- composition of variation into unique and shared components of several sources. For this we can use function varpart of vegan which can handle upto four sources of variation. 
# Function varpart partitions adjusted R2. Unlike ordinary R2 or variance, adjusted R2 is unbiased and its expected value is R2 = 0 for random data (more in lectures). However, it can be easily calculated only for RDA and db-RDA. In principle, CCA can also be used, but calculations are lengthy and complicated, and no software is implemented for them. 

#If the variation is partitioned among groups with the same number of variables (e.g. two soil and two climatic variables), then variation explained by each group is comparable without adjustment. However, if groups contain different numbers of variables, variation explained by not adjusted R2 is not comparable since R2 tends to increase with the number of explanatory variables. Here, use of adjusted R2 is recommended.

#The library vegan offers function varpart, which can partition variation among up to four variables (or groups of variables). Note that varpart is based on redundancy analysis (rda) and uses adjusted R2 to express explained variation. The reason for using only rda is that in R, there is still no function available to calculate adjusted R2 for unimodal ordination methods (like cca). 

###########
# CAUTION
###########
# the below works fine, but I used different matrices by subsetting, which is not what you want persee
# If you compare different datasets, like genetics, morphometrics and environment, this could be of use
# This also works with single vectors as response, which could be of your interest as well!

# make sure "SITES" are depicted as rows

# See detailed documentation:
# vegandocs("partition")

# The first argument gives the dependent data (be it species abundance table, morphological data, or even single vectors (response variable)), and the next two to four arguments give the sources of variation. These can be single variables, matrices or one-sided model formulae like above. 

# Two explanatory matrices -- Hellinger-transform Y
# Formula shortcut "~ ." means: use all variables in 'data'.
mod <- varpart(traits, ~ ., as.matrix(vegdist(traits)), data= data.frame(environment), transfo="hel")
mod

mod <- varpart(traits[,1], ~ ., as.matrix(traits[,2:length(colnames(traits))]), data= data.frame(environment), transfo="hel")
mod
showvarparts(2)
plot(mod)

mod <- varpart(traits[,1], ~ pH + Temp, ~ NH4 + Prot + Carbo, as.matrix(traits[,2:10]), data=environment)


mod <- varpart(traits, ~ pH + Temp, ~ NH4 + Prot + Carbo,data=environment, transfo="hel")


#to interpret the output, read the "Individual fractions", and look at the adjusted R2 (which adds up to 1)

sum(mod$part$indfrac$Adj.R.squared)

mod$part$indfrac

#	how I understand it, is that the residual fraction is your "unexplained" variance
#	come and have a chat if you want more if anything is unclear
library(devtools)
install_github("cran/VennDiagram")
library(VennDiagram)


##### CUSTOM VISUALISATION OF VARPART OBJECTS
# we can use venndiagrams to visualise the data from varpart objects
# for simple varparts (2 indparts), we can also use stacked parplots

install.packages('venneuler')
library(venneuler)


##### FOR 2 sections
#retrieve adjustedRsquared values form the varpart object
Rsq = mod$part$indfrac$Adj.R.squared

#if you want to convert the adj.R.squared values to absolute values (not recommended)
#Rsq = abs(Rsq)

#If you want to change the negative values to 0 (recommended)
Rsq[Rsq<0]=0

size2venn = venneuler(c(A=Rsq[1],B=Rsq[2],"A&B"=Rsq[2]))
size2venn $labels <- c(paste('A',round(Rsq[1],digits=2),sep='\n'),paste('B',round(Rsq[2],digits=2),sep='\n'))
plot(size2venn)
plot(size2venn, col=c("darkred", "forestgreen"),main=paste('unexplained:',round(Rsq[length(Rsq)],digits=3)))

#barcharts for 2 sectors 
df.varpart = data.frame(fraction=c('A','B','A&B','unexplained'),'value'=Rsq,'x'=as.factor(rep(1,length(Rsq))))
levels(df.varpart$fraction) = rev(levels(df.varpart$fraction))
ggplot() + geom_bar(aes(y = value,x=x, fill = fraction), data = df.varpart,
                           stat="identity")+scale_fill_manual(values=c('#D6D6D6','#C6DCBE','#A3A58A','#D0B1B3'))+coord_flip()+ theme_classic()




##### FOR 3 sections
#retrieve adjustedRsquared values form the varpart object
Rsq =mod$part$indfract$Adj.R.square

#if you want to convert the adj.R.squared values to absolute values (not recommended)
#Rsq = abs(Rsq)

#If you want to change the negative values to 0 (recommended)
Rsq[Rsq<0]=0

size3venn <- venneuler(c(A=Rsq[1],B=Rsq[2],C=Rsq[3],"A&B"=Rsq[4], 
                       "A&C"=Rsq[6],"B&C"=Rsq[5],"A&B&C"=Rsq[7]))
size3venn $labels <- c(paste('A',round(Rsq[1],digits=2),sep='\n'),paste('B',round(Rsq[2],digits=2),sep='\n'),paste('C',round(Rsq[3],digits=2),sep='\n'))
plot(size3venn)

df.varpart = data.frame(fraction=c('A','B','C','A&B','B&C','A&C',"A&B&C",'unexplained'),'value'=Rsq,'x'=rep(1,length(Rsq)))
levels(df.varpart$fraction) = rev(levels(df.varpart$fraction))
ggplot() + geom_bar(aes(y = value,x=x, fill = fraction), data = df.varpart,
                           stat="identity")+coord_flip()+ theme_classic()




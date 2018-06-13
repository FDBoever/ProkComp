
# IJSEM table downlaod from 
# https://figshare.com/articles/International_Journal_of_Systematic_and_Evolutionary_Microbiology_IJSEM_phenotypic_database/4272392
# note, this table is rather a mess, and needs extensive cleaning up


#Load the packages used
#install.packages("devtools",dependencies=TRUE)
#library(devtools)
#install_github("vqv/ggbiplot")
library(ggplot2)
library(ggbiplot)

#Reads in the file
RAW = read.table("~/Downloads/4272392/IJSEM_pheno_db_v1.0.txt",header=TRUE,sep="\t")

# table cleanup process
selectedColumns=c(5,7,8,15,17,19)
corrected =c()
for(i in selectedColumns){
	RAW[,i] = gsub("Jan", "1", RAW[,i])
	RAW[,i] = gsub("Feb", "2", RAW[,i])
	RAW[,i] = gsub("Mar", "3", RAW[,i])
	RAW[,i] = gsub("Apr", "4", RAW[,i])
	RAW[,i] = gsub("May", "5", RAW[,i])
	RAW[,i] = gsub("Jun", "6", RAW[,i])
	RAW[,i] = gsub("Jul", "7", RAW[,i])
	RAW[,i] = gsub("Aug", "8", RAW[,i])
	RAW[,i] = gsub("Sep", "9", RAW[,i])
	RAW[,i] = gsub("Oct", "10", RAW[,i])
	RAW[,i] = gsub("Nov", "11", RAW[,i])
	RAW[,i] = gsub("Dec", "12", RAW[,i])

	RAW[,i] = gsub("_01_14", "_1", RAW[,i])
	RAW[,i] = gsub("_02_14", "_2", RAW[,i])
	RAW[,i] = gsub("_03_14", "_3", RAW[,i])
	RAW[,i] = gsub("_04_14", "_4", RAW[,i])
	RAW[,i] = gsub("_05_14", "_5", RAW[,i])
	RAW[,i] = gsub("_06_14", "_6", RAW[,i])
	RAW[,i] = gsub("_07_14", "_7", RAW[,i])
	RAW[,i] = gsub("_08_14", "_8", RAW[,i])
	RAW[,i] = gsub("_09_14", "_9", RAW[,i])
	RAW[,i] = gsub("_10_14", "_10", RAW[,i])
	RAW[,i] = gsub("_11_14", "_11", RAW[,i])
	RAW[,i] = gsub("_12_14", "_12", RAW[,i])
	
	RAW[,i] = gsub("Not_indicated", "", RAW[,i])
	RAW[,i] = gsub("not_indicated", "", RAW[,i])
	RAW[,i] = gsub("n_a", "N", RAW[,i])

	RAW[,i] = gsub("6.75_6.0_7.5", "6.75_7.5", RAW[,i])
	RAW[,i] = gsub("_and_", "_", RAW[,i])
	RAW[,i] = gsub("_9__not_at_6_or_", "_", RAW[,i])
	RAW[,i] = gsub("20g_L", "", RAW[,i])
	RAW[,i] = gsub("<", "", RAW[,i])
	RAW[,i] = gsub(">", "", RAW[,i])
	RAW[,i] = gsub("\x8c\xb10.4", "", RAW[,i])
	RAW[,i] = gsub("\x8c\xb10.5", "", RAW[,i])
	RAW[,i] = gsub("DQ987877", "", RAW[,i])
	values = c()
	for(i in RAW[,i]){
		if(grepl("_", i)){
			min = substr(i, 1, regexpr('_', i)[1]-1)
			max = substr(i, regexpr('_', i)[1]+1, nchar(i))
			print(paste(i, min,max,sep="  - "))
			avg = mean(c(as.numeric(min),as.numeric(max)))
			print(avg)
			values=c(values,avg)

		}else{
			values =c(values,i)
		}
	}
	corrected = cbind(corrected , values)
	}


corrected = data.frame(corrected)
colnames(corrected) = c("GC","meanL","meanW","optPH","optT","optNaCl")

par(mfrow=c(3,2))
corrected$GC = as.numeric(as.character(corrected$GC ))
corrected$GC[corrected$GC<20]=NA
corrected$GC[corrected$GC>80]=NA
hist(as.numeric(corrected$GC ),breaks=60,xlab="%GC",main="")

corrected$meanL = as.numeric(as.character(corrected$meanL ))
corrected$meanL[corrected$meanL>80]=NA
hist(as.numeric(corrected$meanL ),breaks=60,xlab="mean length",main="")

corrected$meanW = as.numeric(as.character(corrected$meanW))
corrected$meanW[corrected$meanW>10]=NA
hist(as.numeric(corrected$meanW),breaks=60,xlab="mean width",main="")

corrected$optPH = as.numeric(as.character(corrected$optPH))
corrected$optPH[corrected$optPH>14]=NA
hist(as.numeric(corrected$optPH),breaks=60,xlab="optimal pH",main="")

corrected$optT = as.numeric(as.character(corrected$optT))
corrected$optT[corrected$optT>121]=NA
hist(as.numeric(corrected$optT),breaks=60,xlab="optimal growth Temperature",main="")

corrected$optNaCl = as.numeric(as.character(corrected$optNaCl))
corrected$optNaCl[corrected$optNaCl>121]=NA
hist(as.numeric(corrected$optNaCl),breaks=60,xlab="optimal c[NaCl]",main="")

selectedColumns=c(16,18,20)
corrected2 =c()
for(i in selectedColumns){
	RAW[,i] = gsub("Jan", "1", RAW[,i])
	RAW[,i] = gsub("Feb", "2", RAW[,i])
	RAW[,i] = gsub("Mar", "3", RAW[,i])
	RAW[,i] = gsub("Apr", "4", RAW[,i])
	RAW[,i] = gsub("May", "5", RAW[,i])
	RAW[,i] = gsub("Jun", "6", RAW[,i])
	RAW[,i] = gsub("Jul", "7", RAW[,i])
	RAW[,i] = gsub("Aug", "8", RAW[,i])
	RAW[,i] = gsub("Sep", "9", RAW[,i])
	RAW[,i] = gsub("Oct", "10", RAW[,i])
	RAW[,i] = gsub("Nov", "11", RAW[,i])
	RAW[,i] = gsub("Dec", "12", RAW[,i])

	RAW[,i] = gsub("_01_14", "_1", RAW[,i])
	RAW[,i] = gsub("_02_14", "_2", RAW[,i])
	RAW[,i] = gsub("_03_14", "_3", RAW[,i])
	RAW[,i] = gsub("_04_14", "_4", RAW[,i])
	RAW[,i] = gsub("_05_14", "_5", RAW[,i])
	RAW[,i] = gsub("_06_14", "_6", RAW[,i])
	RAW[,i] = gsub("_07_14", "_7", RAW[,i])
	RAW[,i] = gsub("_08_14", "_8", RAW[,i])
	RAW[,i] = gsub("_09_14", "_9", RAW[,i])
	RAW[,i] = gsub("_10_14", "_10", RAW[,i])
	RAW[,i] = gsub("_11_14", "_11", RAW[,i])
	RAW[,i] = gsub("_12_14", "_12", RAW[,i])
	
	RAW[,i] = gsub("Not_indicated", "", RAW[,i])
	RAW[,i] = gsub("not_indicated", "", RAW[,i])
	RAW[,i] = gsub("n_a", "N", RAW[,i])

	RAW[,i] = gsub("6.75_6.0_7.5", "6.75_7.5", RAW[,i])
	RAW[,i] = gsub("_to_", "_", RAW[,i])
	RAW[,i] = gsub("_and_", "_", RAW[,i])
	RAW[,i] = gsub("g_L", "", RAW[,i])
	RAW[,i] = gsub("<", "", RAW[,i])
	RAW[,i] = gsub(">", "", RAW[,i])
	RAW[,i] = gsub("_\x8a\xe637___\x8a\xe645", "", RAW[,i])
	RAW[,i] = gsub("___", "_", RAW[,i])
	RAW[,i] = gsub("__", "_", RAW[,i])
	RAW[,i] = gsub("_2014", "", RAW[,i])
	RAW[,i] = gsub("_1937", "", RAW[,i])
	RAW[,i] = gsub("_0_", "", RAW[,i])

	mins = c()
	maxs = c()
	for(i in RAW[,i]){
		print(i)
		if(grepl("_", i)){
			min = substr(i, 1, regexpr('_', i)[1]-1)
			max = substr(i, regexpr('_', i)[1]+1, nchar(i))
			print(paste(i, min,max,sep="  - "))
			mins = c(mins,min)
			maxs = c(maxs,max)
		}else{
			mins = c(mins,"")
			maxs = c(maxs,"")
		}
	}
	corrected2 = cbind(corrected2 , mins,maxs)
	}
	
corrected2 = data.frame(corrected2)
names(corrected2)=c("minPH","maxPH","minT","maxT","minNaCl","maxNaCl")
corrected2$minPH = as.numeric(as.character(corrected2$minPH))
corrected2$maxPH = as.numeric(as.character(corrected2$maxPH))
corrected2$minT = as.numeric(as.character(corrected2$minT))
corrected2$maxT = as.numeric(as.character(corrected2$maxT))
corrected2$minNaCl = as.numeric(as.character(corrected2$minNaCl))
corrected2$maxNaCl = as.numeric(as.character(corrected2$maxNaCl))

cIJSEM = data.frame(cbind(RAW, corrected, corrected2))




########################################################
#
#
#      PLOTS for phylum bacteria 
#
#
########################################################
#FINAL CORRECTED DATATABLE FOR ALL BACTERIA
cIJSEM = data.frame(cbind(RAW, corrected))

p1 = ggplot(cIJSEM, aes(x=optT, y= GC,colour = oxygen_preference)) + geom_point()+theme_classic()
p2 = ggplot(cIJSEM, aes(x= optPH, y= GC)) + geom_point()+theme_classic()
p3 = ggplot(cIJSEM, aes(x= log(optNaCl), y= GC )) + geom_point()+theme_classic()
p4 = ggplot(cIJSEM, aes(x= log(meanL), y= log(meanW))) + geom_point()+theme_classic()
p5 = ggplot(cIJSEM, aes(x= log(meanW), y= GC )) + geom_point()+theme_classic()
p6= ggplot(cIJSEM, aes(x= log(meanW), y= GC)) + geom_point()+theme_classic()

multiplot(p1,p2,p3,p4,p5,p6,cols=2)

p1 = ggplot(cIJSEM, aes(x=optT, y= GC,colour = Habitat)) + geom_point()+theme_classic()
p2 = ggplot(cIJSEM, aes(x= optPH, y= GC,colour = Habitat)) + geom_point()+theme_classic()
p3 = ggplot(cIJSEM, aes(x= log(optNaCl), y= GC ,colour = Habitat)) + geom_point()+theme_classic()
p4 = ggplot(cIJSEM, aes(x= log(meanL), y= log(meanW),colour = Habitat)) + geom_point()+theme_classic()
p5 = ggplot(cIJSEM, aes(x= log(meanW), y= GC ,colour = Habitat)) + geom_point()+theme_classic()
p6= ggplot(cIJSEM, aes(x= log(meanW), y= GC,colour = Habitat)) + geom_point()+theme_classic()
multiplot(p1,p2,p3,p4,p5,p6,cols=2)

p7 = ggplot(cIJSEM, aes(x= oxygen_preference, y= GC)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()
p8 = ggplot(cIJSEM, aes(x= oxygen_preference, y= optT)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5)) +coord_flip()

p9 = ggplot(cIJSEM, aes(x= oxygen_preference, y= optPH)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5)) +coord_flip()

p10 = ggplot(cIJSEM, aes(x= oxygen_preference, y= optNaCl)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()
p11 =ggplot(cIJSEM, aes(x= oxygen_preference, y= meanW)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()
p12 =ggplot(cIJSEM, aes(x= oxygen_preference, y= meanL)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()

multiplot(p7,p8,p9,p10,p11,p12,cols=2)

p7 = ggplot(cIJSEM, aes(x= Habitat, y= GC)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()
p8 = ggplot(cIJSEM, aes(x= Habitat, y= optT)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5)) +coord_flip()
p9 = ggplot(cIJSEM, aes(x= Habitat, y= optPH)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5)) +coord_flip()
p10 = ggplot(cIJSEM, aes(x= Habitat, y= optNaCl)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()
p11 =ggplot(cIJSEM, aes(x= Habitat, y= meanW)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()
p12 =ggplot(cIJSEM, aes(x= Habitat, y= meanL)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()

multiplot(p7,p8,p9,p10,p11,p12,cols=2)


########################################################
########################################################
#
#
#      MAKE A TAXONOMIC SELECTION AND ADD TAXONOMIC INFORMATION
#
#
########################################################
########################################################
# Classificatin information was obtained from bacterio.met (IJSEM 2017 update)
# note, could extract information for the complete set 
BactNet = readLines("~/DATA/MarinobacterGenomics/miscl/Bacterio_net_classifciation 2.txt")

cphylum = c()
cclass=c()
cOrder = c()
cFamily = c()
cgenus = c()

for(line in BactNet){
	list = strsplit(line, split=" ")
	if(list[[1]][1]=="Phylum"){
		phylum = list[[1]][2]
		}else if (list[[1]][1]=="Class"){
			class = list[[1]][2]
		}else if (list[[1]][1]=="Order"){
			order = list[[1]][2]
		}else if (list[[1]][1]=="Family"){
			family = list[[1]][2]
		}else{
			genus=list[[1]]
			genus = data.frame(genus)
			print(genus[1])
			cphylum = c(cphylum,rep(phylum,nrow(genus)))
			cclass = c(cclass,rep(class,nrow(genus)))
			cOrder = c(cOrder,rep(order,nrow(genus)))
			cFamily = c(cFamily,rep(family,nrow(genus)))
			cgenus = c(cgenus, as.character(genus$genus))
	}
}
TaxTable = cbind(Phylum = cphylum,Class= cclass,Order= cOrder,Family= cFamily,Genus_name= cgenus)
TaxTable=data.frame(TaxTable)

cIJSEM_tax <- join(cIJSEM, TaxTable, by = "Genus_name")

ggplot(cIJSEM_tax, aes(x=reorder(Class, GC, FUN=median), y= GC)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+ theme(legend.position="none")+facet_wrap(~ Phylum,scales = "free_x")
ggplot(cIJSEM_tax, aes(x=reorder(Class, maxT-minT, FUN=median), y= maxT-minT,fill=Phylum)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+ theme(legend.position="none")


pcaSET = cIJSEM_tax[,c('Genus_name','Family','Order','Class','Phylum','Habitat','GC','optT','optPH','optNaCl','meanW','meanL','minT','maxT','minPH','maxPH','minNaCl','maxNaCl')]
pcaSET = pcaSET[complete.cases(pcaSET), ]

pcaSET = cIJSEM_tax[,c('Genus_name','Family','Order','Class','Phylum','GC','optT','optPH','optNaCl','minT','maxT','minPH','maxPH','minNaCl','maxNaCl')]
pcaSET = pcaSET[complete.cases(pcaSET), ]


levels(pcaSET $Family) <- c(levels(pcaSET$Family), "Marinobacter_group","Marinobacterium_group")
  levels(pcaSET $Order) <- c(levels(pcaSET$Order), "Marinobacter_group","Marinobacterium_group")
pcaSET$Family[pcaSET$Genus_name == "Marinobacter"] = as.factor("Marinobacter_group")
pcaSET$Order[pcaSET$Genus_name == "Marinobacter"] = "Marinobacter_group"

IJSEM.pca <- prcomp( pcaSET[,7:ncol(pcaSET)], scale. = TRUE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,2], ellipse = TRUE)+ scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')


pcaSET = subset(pcaSET,Phylum=="Acidobacteria")
pcaSET = subset(pcaSET,Class=="Gammaproteobacteria")
pcaSET = subset(pcaSET,Order=="Oceanospirillales")

ggplot(pcaSET, aes(x=reorder(Genus_name, maxT-minT, FUN=median), y= maxT-minT,fill= Genus_name)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+ theme(legend.position="none")
ggplot(pcaSET, aes(x=reorder(Habitat, maxT-minT, FUN=median), y= maxT-minT,fill= Habitat)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+ theme(legend.position="none")
ggplot(pcaSET, aes(x=reorder(Habitat, maxNaCl-minNaCl, FUN=median), y=maxNaCl-minNaCl,fill= Habitat)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+ theme(legend.position="none")
ggplot(pcaSET, aes(x=reorder(Habitat, maxPH-minPH, FUN=median), y=maxPH-minPH,fill= Habitat)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+ theme(legend.position="none")


IJSEM.pca <- prcomp( pcaSET[,6:ncol(pcaSET)], scale. = TRUE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,3], ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

IJSEM.pca <- prcomp( pcaSET[,6:ncol(pcaSET)], scale. = TRUE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,3]) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

IJSEM.pca <- prcomp( pcaSET[,6:ncol(pcaSET)], scale. = TRUE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,2]) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

IJSEM.pca <- prcomp( pcaSET[,6:ncol(pcaSET)], scale. = TRUE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,2]) + scale_color_manual(values=getPalette(length(unique(pcaSET[,2]))))+
  theme(legend.direction = 'horizontal', legend.position = 'top')+theme_classic()

IJSEM.pca <- prcomp( pcaSET[,6:ncol(pcaSET)], scale. = TRUE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,3],ellipse = TRUE) + scale_color_manual(values=getPalette(length(unique(pcaSET[,3]))))+
  theme(legend.direction = 'horizontal', legend.position = 'top')+theme_classic()

pcaSET = subset(pcaSET,Order %in% c('Marinobacter_group','Alteromonadales','Oceanospirillales'))

IJSEM.pca <- prcomp( pcaSET[,6:ncol(pcaSET)], scale. = TRUE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,3],ellipse = TRUE) + scale_color_manual(values=getPalette(length(unique(pcaSET[,3]))))+
  theme(legend.direction = 'horizontal', legend.position = 'top')+theme_classic()

IJSEM.pca <- prcomp( pcaSET[,6:ncol(pcaSET)], scale. = FALSE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,2],ellipse = TRUE) + scale_color_manual(values=getPalette(length(unique(pcaSET[,2]))))+
  theme(legend.direction = 'horizontal', legend.position = 'top')+theme_classic()

IJSEM.pca <- prcomp( pcaSET[,6:ncol(pcaSET)], scale. = FALSE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,3],ellipse = TRUE) + scale_color_manual(values=getPalette(length(unique(pcaSET[,3]))))+
  theme(legend.direction = 'horizontal', legend.position = 'top')+theme_classic()

IJSEM.pca <- prcomp( pcaSET[,6:ncol(pcaSET)], scale. = TRUE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,3],ellipse = TRUE) + scale_color_manual(values=getPalette(length(unique(pcaSET[,3]))))+
  theme(legend.direction = 'horizontal', legend.position = 'top')+theme_classic()

pcaSET = subset(pcaSET,Order %in% c('Marinobacter_group','Oceanospirillales'))
IJSEM.pca <- prcomp( pcaSET[,6:ncol(pcaSET)], scale. = TRUE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,2],ellipse = TRUE) + scale_color_manual(values=getPalette(length(unique(pcaSET[,2]))))+
  theme(legend.direction = 'horizontal', legend.position = 'top')+theme_classic()


 gAlter = c("Aestuariibacter","Agaribacter","Agarivorans","Aliagarivorans","Aliiglaciecola","Alishewanella","Alteromonas","Bowmanella","Catenovulum","Glaciecola","Marinobacter","Marinobacterium","Melitea","Salinimonas","Tamilnaduibacter","Celerinatantimonas","Colwellia","Thalassomonas","Thalassotalea","Ferrimonas","Paraferrimonas","Idiomarina","Pseudidiomarina","Moritella","Paramoritella","Algicola","Pseudoalteromonas","Psychrosphaera","Psychromonas","Psychrobium","Shewanella")	
 
 fAlter = c(rep("Alteromonadaceae",15),c("Celerinatantimonadaceae","Colwelliaceae","Colwelliaceae","Colwelliaceae","Ferrimonadaceae","Ferrimonadaceae","Idiomarinaceae","Idiomarinaceae","Moritellaceae","Moritellaceae","Pseudoalteromonadaceae","Pseudoalteromonadaceae","Pseudoalteromonadaceae","Psychromonadaceae","Shewanellaceae","Shewanellaceae"))
 
 fOceano = c("Saccharospirillaceae","Saccharospirillaceae","Saccharospirillaceae","Alcanivoracaceae","Alcanivoracaceae","Alcanivoracaceae","Alcanivoracaceae","Hahellaceae","Hahellaceae","Hahellaceae","Hahellaceae","Hahellaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Halomonadaceae","Litoricolaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oceanospirillaceae","Oleiphilaceae","unassigned","unassigned","unassigned")

gOceano = c("Gynuella","Saccharospirillum","Salinispirillum","Alcanivorax","Fundibacter","Kangiella","Pleionea","Endozoicomonas","Hahella","Halospina","Kistimonas","Zooshikella","Aidingimonas","Carnimonas","Chromohalobacter","Cobetia","Deleya","Halomonas","Halotalea","Halovibrio","Kushneria","Larsenimonas","Modicisalibacter","Salinicola","Terasakiispira","Volcaniella","Zymobacter","Litoricola","Amphritea","Balneatrix","Bermanella","Corallomonas","Litoribrevibacter","Marinomonas","Marinospirillum","Neptuniibacter","Neptunomonas","Nitrincola","Oceaniserpentilla","Oceanobacter","Oceanospirillum","Oleibacter","Oleispira","Pseudospirillum","Reinekea","Thalassolituus","Oleiphilus","Motiliproteus","Salicola","Spongiispira")
 
 
cIJSEM_select = subset(cIJSEM, Genus_name %in% c(gAlter,gOceano))

TaxTable = cbind(Genus_name = c(gAlter, gOceano),family=c(fAlter, fOceano),order=c(rep("Alteromonadales",length(fAlter)), rep("Oceanospirillales",length(fOceano))))

TaxTable=data.frame(TaxTable)

cIJSEM_select <- join(cIJSEM_select, TaxTable, by = "Genus_name")


q1 = ggplot(cIJSEM_select, aes(x=reorder(Genus_name, GC, FUN=median), y= GC,fill=order)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+facet_wrap( ~ order,scales = "free_x")+ theme(legend.position="none")

q2 =ggplot(cIJSEM_select, aes(x=reorder(Genus_name, optT, FUN=median), y= optT,fill=order)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+facet_grid( ~ order,scales = "free_x")+ theme(legend.position="none")
q3 =ggplot(cIJSEM_select, aes(x=reorder(Genus_name, optPH, FUN=median), y= optPH,fill=order)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+facet_grid( ~ order,scales = "free_x")+ theme(legend.position="none")
q4 =ggplot(cIJSEM_select, aes(x=reorder(Genus_name, optNaCl, FUN=median), y= optNaCl,fill=order)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+facet_grid( ~ order,scales = "free_x")+ theme(legend.position="none")
q5 =ggplot(cIJSEM_select, aes(x=reorder(Genus_name, meanW, FUN=median), y= meanW,fill=order)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+facet_grid( ~ order,scales = "free_x")+ theme(legend.position="none")
q6 =ggplot(cIJSEM_select, aes(x=reorder(Genus_name, meanL, FUN=median), y= meanL,fill=order)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+facet_grid( ~ order,scales = "free_x")+ theme(legend.position="none")

multiplot(q1,q2,q3,q4,q5,q6,cols=3)

ggplot(cIJSEM_select, aes(x=optNaCl, y= GC,fill=order)) + geom_point(size=2,shape=21, aes(fill = factor(order))) + scale_fill_manual(values=c("chocolate3", "cyan4"))+theme_classic()

ggplot(cIJSEM_select, aes(x=reorder(Genus_name, GC), y= GC, fill= order)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()
ggplot(cIJSEM_select, aes(x=reorder(Genus_name, optNaCl), y= optNaCl, fill= order)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()

ggplot(cIJSEM_select, aes(x=reorder(Genus_name, order), y= meanW, fill= order)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()

ggplot(cIJSEM_select, aes(x=reorder(Genus_name, GC), y= GC)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))
ggplot(cIJSEM_select, aes(x=reorder(Genus_name, optT), y= optT)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))
ggplot(cIJSEM_select, aes(x=reorder(Genus_name, optPH), y= optPH)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))
ggplot(cIJSEM_select, aes(x=reorder(Genus_name, optNaCl), y= optNaCl)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))
ggplot(cIJSEM_select, aes(x=reorder(Genus_name, optNaCl), y= optNaCl)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))


ggplot(cIJSEM_select, aes(x=reorder(Genus_name, GC), y= GC, fill= orde)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))
ggplot(cIJSEM_select, aes(x=reorder(family, GC), y= GC, fill= family)) + geom_boxplot() + xlab(" ")+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5))+coord_flip()

########################################################
# Principal Component Analysis and BIPLOTS 
########################################################


pcaSET = cIJSEM_select[,c('Genus_name','family','order','GC','optT','optPH','optNaCl','meanW','meanL')]
pcaSET = pcaSET[complete.cases(pcaSET), ]

levels(pcaSET $family) <- c(levels(pcaSET$family), "Marinobacter_group","Marinobacterium_group")
  levels(pcaSET $order) <- c(levels(pcaSET$order), "Marinobacter_group","Marinobacterium_group")
pcaSET$family[pcaSET$Genus_name == "Marinobacter"] = as.factor("Marinobacter_group")
pcaSET$order[pcaSET$Genus_name == "Marinobacter"] = "Marinobacter_group"


IJSEM.pca <- prcomp( pcaSET[,4:ncol(pcaSET)], scale. = TRUE)
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,3], ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,2], ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
  
ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,2], ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

#without morphology
IJSEM.pca <- prcomp( pcaSET[,4:(ncol(pcaSET)-2)], scale. = TRUE)
b1 = ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,3], ellipse = TRUE, circle = FALSE) +
  scale_color_discrete(name = '') #+
  theme(legend.direction = 'horizontal', legend.position = 'top')

b2 = ggbiplot(IJSEM.pca, obs.scale = 1, var.scale = 1,
  groups = pcaSET[,2], ellipse = TRUE, circle = FALSE) +
  scale_color_discrete(name = '') #+
  theme(legend.direction = 'horizontal', legend.position = 'top')
  
 
pcaSET2 = subset(pcaSET, order  %in% c('Alteromonadales','Marinobacter_group'))

IJSEM.pca2 <- prcomp( pcaSET2[,4:ncol(pcaSET2)], scale. = TRUE)
ggbiplot(IJSEM.pca2, obs.scale = 1, var.scale = 1,
  groups = pcaSET2[,2], ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

# without morphology
IJSEM.pca2 <- prcomp( pcaSET2[,4:(ncol(pcaSET2)-2)], scale. = TRUE)
b3=ggbiplot(IJSEM.pca2, obs.scale = 1, var.scale = 1,
  groups = pcaSET2[,2], ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '')# +
  theme(legend.direction = 'horizontal', legend.position = 'top')

pcaSET3 = subset(pcaSET, order  %in% c('Oceanospirillales','Marinobacter_group'))

IJSEM.pca3 <- prcomp( pcaSET3[,4:ncol(pcaSET3)], scale. = TRUE)
ggbiplot(IJSEM.pca3, obs.scale = 1, var.scale = 1,
  groups = pcaSET3[,2], ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

# without morphology
IJSEM.pca3 <- prcomp( pcaSET3[,4:(ncol(pcaSET3)-2)], scale. = TRUE)
b4=ggbiplot(IJSEM.pca3, obs.scale = 1, var.scale = 1,
  groups = pcaSET3[,2], ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') #+
  theme(legend.direction = 'horizontal', legend.position = 'top')
  
multiplot(b1,b2,b3,b4,cols=2)

hist(pcaSET3$optT),breaks=60,xlab="optimal pH",main="")
hist(pcaSET3$optPH,breaks=60,xlab="optimal pH",main="")
hist(pcaSET3$optNaCl,breaks=60,xlab="optimal pH",main="")
hist(pcaSET3$GC,breaks=60,xlab="optimal pH",main="")

h1 = ggplot(pcaSET, aes(x=GC, color=order,fill=order)) +
  geom_histogram( alpha=0.5, position="identity", binwidth=1)
h2 = ggplot(pcaSET, aes(x=optT, color=order,fill=order)) +
  geom_histogram( alpha=0.5, position="identity",binwidth=1)
  h3 = ggplot(pcaSET, aes(x=optPH, color=order,fill=order)) +
  geom_histogram( alpha=0.5, position="identity",binwidth=0.25)
  h4 = ggplot(pcaSET, aes(x=optNaCl, color=order,fill=order)) +
  geom_histogram( alpha=0.5, position="identity")
  
  h5 = ggplot(pcaSET, aes(x=meanW, color=order,fill=order)) +
  geom_histogram( alpha=0.5, position="identity",binwidth=0.25)
  h6 = ggplot(pcaSET, aes(x= meanL, color=order,fill=order)) +
  geom_histogram( alpha=0.5, position="identity")
multiplot(h1,h2,h3,h4,h5,h6,cols=2)

geom_density
h1 = ggplot(pcaSET, aes(x=GC, color=order,fill=order)) +
  geom_density( alpha=0.5, position="identity", binwidth=1)
h2 = ggplot(pcaSET, aes(x=optT, color=order,fill=order)) +
  geom_density( alpha=0.5, position="identity",binwidth=1)
  h3 = ggplot(pcaSET, aes(x=optPH, color=order,fill=order)) +
  geom_density( alpha=0.5, position="identity",binwidth=0.25)
  h4 = ggplot(pcaSET, aes(x=optNaCl, color=order,fill=order)) +
  geom_density( alpha=0.5, position="identity")
  
  h5 = ggplot(pcaSET, aes(x=meanW, color=order,fill=order)) +
  geom_density( alpha=0.5, position="identity",binwidth=0.25)
  h6 = ggplot(pcaSET, aes(x= meanL, color=order,fill=order)) +
  geom_density( alpha=0.5, position="identity")
multiplot(h1,h2,h3,h4,h5,h6,cols=2)

hist(pcaSET3$optT,break=60)


########################################################
# WRITE THE DATATABLES AS FILES
########################################################


write.csv(cIJSEM,"~/IJSEM_phenotype_FDB.txt")
write.csv(pcaSET,"~/IJSEM_Altero_Oceano_FDB.txt")
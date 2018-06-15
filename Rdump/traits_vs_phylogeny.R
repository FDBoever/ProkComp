library(devtools)
devtools::install_github("fkeck/phylosignal")
devtools::install_github("fmichonneau/phylobase")
install.packages("phylosignal",dependencies=TRUE)

install.packages(c('adephylo','phylobase','picante'),dependencies=TRUE)
library(phylosignal)
library(adephylo)
library(ape)
library(phylobase)
library(phytools)
library(picante)
#______________________________________________________________________________________________________________


RAxMLANVIO = read.tree("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/ANVIO_RAxML/RAxML_bipartitions.autosubst100")
RAxMLANVIORooted =  midpoint.root(RAxMLANVIO)
tre <- ladderize(RAxMLANVIORooted, right = FALSE)

tre$node.label <- NULL # Erase the bootstrap values from the phylo object

CheckM = read.table("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/QC/qa2Marinobacter.txt",header=TRUE)
CheckM2 = subset(CheckM, Completeness > 98)
CheckM2 = CheckM2[order(CheckM2$Bin_Id),]
CheckM2 = CheckM2[CheckM2$Bin_Id != "Marinobacter_dagiaonensis_CGMCC_1_9167",]

is_tip <- tre $edge[,2] <= length(tre $tip.label)
ordered_tips <- tre $edge[is_tip, 2]
tre$tip.label[ordered_tips]

#order dataframe according to phylogeny
CheckM2 = CheckM2[match(tre$tip.label[ordered_tips], CheckM2$Bin_Id),]
CheckM2 <- CheckM2[seq(dim(CheckM2)[1],1),]


dat <- list()
dat$GC <- CheckM2$GC
dat$random <- rnorm(60, sd = 5)
dat$bm <- rTraitCont(tre)
dat <- as.data.frame(dat)

p4d <- phylo4d(tre, dat)
p5d = data.frame(p4d)
barplot.phylo4d(p4d, tree.type = "phylo", tree.ladderize = TRUE)
dotplot(p4d,tree.ladderize = TRUE)
barplot(p4d, trait = c("bm", "GC"))

barplot(p4d, tree.type = "fan", tip.cex = 0.6, tree.open.angle = 180, trait.cex = 0.6)

mat.col <- ifelse(tdata(p4d, "tip") < 0, "red", "grey35")
barplot(p4d, tree.ladderize = TRUE,bar.col = mat.col)

tip.col <- rep(1, nTips(p4d))
tip.col[(mat.col[, 2] == "red") | (mat.col[, 3] == "red")] <- 2
barplot(p4d, tree.ladderize = TRUE,trait.bg.col = c("#F6CED8", "#CED8F6", "#CEF6CE"),
        bar.col = mat.col, tip.col = tip.col, trait.font = c(1, 2, 2))


# Assessing the behavior of these methods with this phylogeny along a Brownian-Motion influence gradient
phylosim <- phyloSim(tree = tre, method = "all", nsim = 100, reps = 99)
plot(phylosim, stacked.methods = FALSE, quantiles = c(0.05, 0.95))
plot.phylosim(phylosim, what = "pval", stacked.methods = TRUE)

#Assessing the signal depth with correlograms

mass.crlg <- phyloCorrelogram(p4d, trait = "mass")
random.crlg <- phyloCorrelogram(p4d, trait = "random")
bm.crlg <- phyloCorrelogram(p4d, trait = "bm")

plot(mass.crlg)
plot(random.crlg)
plot(bm.crlg)

#Locating the signal with LIPA
carni.lipa <- lipaMoran(p4d[1:60,])
carni.lipa.p4d <- lipaMoran(p4d, as.p4d = TRUE)
barplot.phylo4d(p4d, bar.col=(carni.lipa$p.value < 0.05) + 1, center = FALSE , scale = FALSE)
barplot.phylo4d(carni.lipa.p4d, bar.col = (carni.lipa$p.value < 0.05) + 1, center = FALSE, scale = FALSE)

#______________________________________________________________________________________________________________

#Phylogenetic signal means that closely related species have similar traits. This violates the assumption of independence of data points that is inherent in many methods including correlation and regression (Felsenstein 1985). We can account for non-independence due to phylogenetic signal using methods including phylogenetically independent contrasts and phylogenetic generalised least squares (pGLS).
#Generalised least squares methods work just like an ANOVA or linear model - we can test for relationships between categorical or continuous values, optionally taking phylogenetic relatedness into account.
#Let's test for a relationship between specific root length (SRL) and root tissue density, taking phylogenetic relationships among species into account.
#______________________________________________________________________________________________________________


phy <- read.tree("~/DATA/MarinobacterGenomics/FINAL_PROKCOMP/ANVIO_RAxML/RAxML_bipartitions.autosubst100")
class(phy)
phy

phy = tre


rownames(CheckM2) = CheckM2$Bin_Id
CheckM2 = CheckM2[,2:ncol(CheckM2)]
CheckM2 = CheckM2[,c('Genome_size','GC','GC_sd','coding_density','predicted_genes')]


# check for mismatches/missing species
combined <- match.phylo.data(phy, CheckM2)
# the resulting object is a list with $phy and $data elements.  replace our
# original data with the sorted/matched data
phy <- combined$phy
traits <- combined$data

apply(traits, 2, Kcalc, phy)
multiPhylosignal(traits, multi2di(phy))


plot(phy, direction = "up", show.tip.label = FALSE, show.node.label = TRUE, 
    cex = 0.7)
# Plot leaf area on the phylogeny. cex argument scales symbol size by trait
# value.
tiplabels(pch = 19, col = "black", cex = 3 * traits[, "GC"]))




# GLS of coding_density as a function of GC - non-phylogenetic model
root.gls <- gls(coding_density ~ GC, data = traits)
anova(root.gls)

root.pgls <- gls(coding_density ~ GC, correlation = corBrownian(value = 1, 
    phy), data = traits)
anova(root.pgls)


plot(coding_density ~ GC, data = traits, xlab = "GC ", 
    ylab = "coding_density")
# add model fit lines - coef is the model fit coefficients, lwd increases
# line width
abline(coef(root.gls), lwd = 2, col = "black")
abline(coef(root.pgls), lwd = 2, col = "red")
legend("bottomleft", legend = c("GLS fit", "Phylogenetic GLS fit"), lwd = 2, 
    col = c("black", "red"))





# GLS of predicted_genes as a function of Genome_size - non-phylogenetic model
root.gls <- gls(predicted_genes ~ Genome_size, data = traits)
anova(root.gls)

# Phylogenetic GLS - adds effect of phylogeny to the model
root.pgls <- gls(predicted_genes ~ Genome_size, correlation = corBrownian(value = 1, 
    phy), data = traits)
anova(root.pgls)


plot(predicted_genes ~ Genome_size, data = traits, xlab = "Genome_size ", 
    ylab = "predicted_genes",pch=21)
# add model fit lines - coef is the model fit coefficients, lwd increases
# line width
abline(coef(root.gls), lwd = 2, col = "black")
abline(coef(root.pgls), lwd = 2, col = "red")
legend("topleft", legend = c("GLS fit", "Phylogenetic GLS fit"), lwd = 2, 
    col = c("black", "red"))
    

    
# GLS of predicted_genes as a function of Genome_size - non-phylogenetic model
root.gls <- gls(GC ~ Genome_size, data = traits)
anova(root.gls)

# Phylogenetic GLS - adds effect of phylogeny to the model
root.pgls <- gls(GC  ~ Genome_size, correlation = corBrownian(value = 1, 
    phy), data = traits)
anova(root.pgls)


plot(GC  ~ Genome_size, data = traits, xlab = "Genome_size ", 
    ylab = "GC",pch=21)
# add model fit lines - coef is the model fit coefficients, lwd increases
# line width
abline(coef(root.gls), lwd = 2, col = "black")
abline(coef(root.pgls), lwd = 2, col = "red")
legend("topleft", legend = c("GLS fit", "Phylogenetic GLS fit"), lwd = 2, 
    col = c("black", "red"))

 #There is a weak relationship between genome-size and GC content. The relationship is borderline significant if we do not take phylogenetic relatedness into account. We see a stronger and significant relationship between GC and genome size after taking phylogenetic relatedness into account.


 install.packages('caper',dependencies=TRUE)
 library(caper)
 
 
 
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
mod1 <- pgls(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird, lambda='ML')
mod2 <- pgls(log(Egg.Mass) ~ log(M.Mass), data=shorebird, lambda='ML', delta='ML')

data(BritishBirds)
BritishBirds <- comparative.data(BritishBirds.tree, BritishBirds.data, binomial)
redPhyloD <- phylo.d(BritishBirds, binvar=Red_list)
print(redPhyloD) 
traits$polar = c(rep(0,10),rep(1,8),rep(0,60-18))
traits$names = rownames(traits)
CtraitTree <- comparative.data(phy, traits,names)
phyloD <- phylo.d(CtraitTree, binvar= polar )

print(phyloD)
plot(phyloD) #Scaling of D values for the observed sum of character change from the distributions of simulated sums of character change under models of random association (blue line) and a Brownian process (red line). The distributions of the simulations are used to test the significance of departures from either model.


clmat <- clade.matrix(tre)
pd.bootstrap(clmat, ntips=59, reps=1000, method='TBL')


CM <- clade.matrix(tre)
ED <- ed.calc(CM)
install.packages(c('pGLS','surface','bayou','OUwie'),dependencies=TRUE)

library(pGLS)
library(geiger)
library(caper)
library(picante)
library(surface)
library(igraph)
library(phytools)
library(bayou)
library(OUwie)

tree <- read.tree(tree_file_name)
tree = tre
plot.phylo(tree)
nodelabels()

data <- read.table(data_file_name, row.names=1)
traits <- data[[1]]
names(traits) <- rownames(data)
stderrs <- data[[2]]
names(stderrs) <- rownames(data)


TRAITS=traits
traits =TRAITS[[2]]
names(traits) <- rownames(TRAITS)

stderrs <- TRAITS[[3]]
names(stderrs) <- rownames(TRAITS)


# # Use the phylosignal function of the picante package to test for phylogenetic signal, measured with the K statistic.
# K ~ 1.0 suggests Brownian motion, K < 1 indicates less resemblance among relatives than expected under Brownian motion,
# and K > 1 indicates more resemblance among relatives than expected under Brownian motion (Blomberg et al., 2003, Evolution, 57:717-745).

phylosignal(traits, tree)

# Use the phylosig function of the phytools package to test for phylogenetic signal, measured with Pagel's lambda.
# Pagel's lambda is 0 <= lambda <= 1, with lambda ~ 0 indicating no phylogenetic signal, and lambda ~ 1 indicating as much
# phylogenetic signal as expected under Brownian motion (Pagel, 1999, Nature, 401:877-884).

phylosig(tree, traits, method="lambda", test=TRUE)
phylo.signal.disc(tree, traits, method="lambda", test=TRUE)

# Compare the fit of standard single-regime models, using geiger.

bm.model <- fitContinuous(tree, traits, model="BM")
white.model <- fitContinuous(tree, traits, model="white")
ou.model <- fitContinuous(tree, traits, model="OU")
eb.model <- fitContinuous(tree, traits, model="EB")
results <- c(bm.model$opt$aicc,white.model$opt$aicc,ou.model$opt$aicc,eb.model$opt$aicc)
names(results) <- c("bm","white_noise","ou","eb")
results


prior <- make.prior(tree, dists=list(dalpha="dlnorm", dsig2="dlnorm", dsb="dsb", dk="cdpois",
	dtheta="dnorm"), param=list(dalpha=list(meanlog=-5, sdlog=2),
	dsig2=list(meanlog=-1, sdlog=2), dk=list(lambda=6, kmax=200),
	dsb=list(bmax=Inf,prob=tree$edge.length), dtheta=list(mean=mean(traits), sd=2)))
	
fit1 <- bayou.mcmc(tree, traits, SE=stderrs, model="OU", prior, ngen=2000000, new.dir=TRUE, plot.freq=400000, ticker.freq=200000)


chain1 <- load.bayou(fit1, save.Rdata=FALSE, cleanup=TRUE)
chain1 <- set.burnin(chain1, 0.3)
out <- summary(chain1)

# Reset to use one plots per figure.
par(mfrow=c(1,1))

# Produce a tree figure for the first bayou mcmc chain.
plotSimmap.mcmc(tree, chain1, burnin=0.3, circle=TRUE, fsize=0.4)

# Produce a phenogram for the first bayou mcmc chain.
phenogram.density(tree, traits, chain=chain1, burnin=0.3, pp.cutoff=0.1)
plot(chain1)
# Run a second bayou mcmc chain.
fit2 <- bayou.mcmc(tree, traits,SE=stderrs, model="OU", prior, ngen=2000000, new.dir=TRUE, plot.freq=NULL, ticker.freq=200000)

# Summarize the second bayou mcmc chain.
chain2 <- load.bayou(fit2, save.Rdata=FALSE, cleanup=FALSE)
chain2 <- set.burnin(chain2, 0.3)

true.pars <- list(alpha=1, sig2=1, beta1=0.75, k=3, ntheta=4, theta=c(2,3,1,8))
shifts <- identifyBranches(tree, true.pars$k)
true.pars$sb <- shifts$sb
true.pars$loc <- shifts$loc
true.pars$t2 <- 2:true.pars$ntheta
dat <- dataSim(true.pars, model = "OU", tree, SE=0)$dat

cache <- bayou:::.prepare.ou.univariate(tree, traits,SE=0, pred=pred)
nn <- names(dat)
dat <- dat + true.pars$beta1*cache$pred
dat <- dat[,1]
names(dat) <-nn

postburn <- round((0.3*length(chain1$gen)),0):length(chain1$gen)
est.pars <- list(alpha=median(chain1$alpha[postburn]), sig2=median(chain1$sig2[postburn]), beta1=median(chain1$beta1[postburn]))
est.pars$sb <- which(sumchain$branch.posteriors[,1] > 0.5)
est.pars$k <- length(est.pars$sb)
est.pars$ntheta <- est.pars$k+1
est.pars$loc <- sumchain$branch.posteriors[est.pars$sb,4]
est.pars$t2 <- 1:length(est.pars$sb)+1
est.pars$theta <- c(median(sapply(chain$theta[postburn], function(x) x[1])), sumchain$branch.posteriors[est.pars$sb,2])




#########

haemulidTrees <- read.nexus(url("http://schmitzlab.info/Trees4dryad.nex"))
haemulidTrees
haemulids <- read.csv(url("http://schmitzlab.info/haemulids.csv"), header=TRUE)
head(haemulids)
#Grab the tree from above
base.tree <- haemulidTreesLadderized[[1]]

#assign row names to the dataframe
rownames(haemulids) <- haemulids$Taxon


fit.ER <- ace(z, tree, model="ER", type="discrete")
fit.ER


#divide all branches by 100
scaled.tree <- tree
scaled.tree$edge.length <- scaled.tree$edge.length/100

#try again...
#calculating marginal likelihoods for ancestral states with an all-rates-equal model (ER)
fit.ER <- ace(z, scaled.tree, model="ER", type="discrete")
fit.ER


par(mar=c(1,0,0,1))

#plot tree
plot(tree, label.offset=1.5, cex=0.5); add.scale.bar()

#add tiplabels
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

#add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:tree$Nnode+Ntip(tree), pie=fit.ER$lik.anc, cex=0.35, piecol=c("red", "blue"))

sim.hab <- make.simmap(scaled.tree, z, model="ER")
cols <- setNames(c("blue","red"),c("reef", "non-reef"))
plotSimmap(sim.hab, cols, fsize=0.7, ftype="i")

sim.hab25 <- make.simmap(scaled.tree, z, model="ER", nsim=25)

par(mfrow=c(5,5))
for (i in 1:25){
plotSimmap(sim.hab25[[i]], cols, ftype="off")
}

#make 100 simulations of trait histories
sim.hab100 <- make.simmap(scaled.tree, z, model="ER", nsim=100)


sim.summary <- describe.simmap(sim.hab100, plot=FALSE)
sim.summary

par(mar=c(1,0,0,1))
plot(tree, label.offset=1.5, cex=0.5); add.scale.bar()

#add tiplabels
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

#add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:tree$Nnode+Ntip(tree), pie=sim.summary$ace, cex=0.35, piecol=c("red", "blue"))









denitrif = c('part','part','no','no','no','part','part','no','no','no','no','no','no','no','no','no','no','no','part','yes','yes','yes','yes','yes','yes','yes','yes','part','part','part','no','part','part','part','part','part','part','part','part','no','no','part','part','part','part','yes','part','part','part','part','yes','yes','part','no','no','no','part','part','part','part')
names(denitrif) = tre$tip.label[ordered_tips]
z=denitrif

mycol<-character(length(z))
mycol[z=="yes"]<-"darkolivegreen2"
mycol[z=="no"]<-"brown2"
mycol[z=="part"]<-"darkgoldenrod1"

plot(tre, label.offset=1.5, cex=0.5)
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

write.tree(tre, "haemulid.tre")
tree <- read.tree("haemulid.tre")

tr2<-tree
tr2$edge.length<-NULL
par(mar=c(1,0,0,1))

plot(tr2, label.offset=1.5, cex=0.5)
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)
fit.ER <- ace(z, tree, model="ER", type="discrete")
fit.ER

scaled.tree <- tree
scaled.tree$edge.length <- scaled.tree$edge.length/100
fit.ER <- ace(z, scaled.tree, model="ER", type="discrete")
fit.ER


#plot tree
plot(tr2, label.offset=1.5, cex=0.5); add.scale.bar()
#add tiplabels
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

#add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:tree$Nnode+Ntip(tree), pie=fit.ER$lik.anc, cex=0.35, piecol=c("brown2",'darkgoldenrod1','darkolivegreen2'))
sim.hab <- make.simmap(scaled.tree, z, model="ER")

cols <- setNames(c("red","yellow",'green'),c("no", "part",'yes'))
plotSimmap(sim.hab, cols, fsize=0.7, ftype="i")
sim.hab25 <- make.simmap(scaled.tree, z, model="ER", nsim=25)

par(mfrow=c(5,5))
 for (i in 1:25){
 plotSimmap(sim.hab25[[i]], cols, ftype="off")
}
sim.hab100 <- make.simmap(scaled.tree, z, model="ER", nsim=100)

sim.summary <- describe.simmap(sim.hab100, plot=FALSE)
sim.summary

par(mfrow=c(1,1))

 #proceed to plotting
par(mar=c(1,0,0,1))
plot(tr2, label.offset=1.5, cex=0.5); add.scale.bar()

#add tiplabels
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

#add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:tree$Nnode+Ntip(tree), pie=sim.summary$ace, cex=0.35, piecol=c("brown2",'darkgoldenrod1','darkolivegreen2'))


############


denitrif = c('yes','yes','no','no','no','yes','yes','no','no','no','no','no','no','no','no','no','no','no','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','no','yes','yes','yes','yes','yes','yes','yes','yes','no','no','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','yes','no','no','no','yes','yes','yes','yes')
names(denitrif) = tre$tip.label[ordered_tips]
z=denitrif

mycol<-character(length(z))
mycol[z=="yes"]<-"darkolivegreen2"
mycol[z=="no"]<-"brown2"


plot(tre, label.offset=1.5, cex=0.5)
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

write.tree(tre, "haemulid.tre")
tree <- read.tree("haemulid.tre")

tr2<-tree
tr2$edge.length<-NULL
par(mar=c(1,0,0,1))

plot(tr2, label.offset=1.5, cex=0.5)
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)
fit.ER <- ace(z, tree, model="ER", type="discrete")
fit.ER

scaled.tree <- tree
scaled.tree$edge.length <- scaled.tree$edge.length/100
fit.ER <- ace(z, scaled.tree, model="ER", type="discrete")
fit.ER


#plot tree
plot(tr2, label.offset=1.5, cex=0.5); add.scale.bar()
#add tiplabels
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.3)

#add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:tree$Nnode+Ntip(tree), pie=fit.ER$lik.anc, cex=0.35, piecol=c("brown2",'darkolivegreen2'))
sim.hab <- make.simmap(scaled.tree, z, model="ER")

cols <- setNames(c("red","yellow",'green'),c("no", "part",'yes'))
plotSimmap(sim.hab, cols, fsize=0.7, ftype="i")
sim.hab25 <- make.simmap(scaled.tree, z, model="ER", nsim=25)

par(mfrow=c(5,5))
 for (i in 1:25){
 plotSimmap(sim.hab25[[i]], cols, ftype="off")
}
sim.hab100 <- make.simmap(scaled.tree, z, model="ER", nsim=100)

sim.summary <- describe.simmap(sim.hab100, plot=FALSE)
sim.summary


 #proceed to plotting
par(mar=c(1,0,0,1))
plot(tr2, label.offset=5, cex=0.5); add.scale.bar()

#add tiplabels
tiplabels(pch=21, col="black", adj=3, bg=mycol, cex=1.3)

#add nodelabels (pie charts representing likelihoods)
nodelabels(node=1:tree$Nnode+Ntip(tree), pie=sim.summary$ace, cex=1, piecol=c("brown2",'darkolivegreen2'))
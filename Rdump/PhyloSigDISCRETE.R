# Phylogenetic signal test for discrete trait
#--------------------------------------------------

#----	RATIONALE		-------------------------
#For discrete variable (isolation habitat), we compare 3 models (ER, SYM, ARD) to fit the trait data. SYN was accepted as the best model, we then used the same function with transform=’lamda’ to get the likelihood based on the original tree and compared it with the likelihood based on lambda 0 tree, we used a likelihood ratio test to assess wheter the trait shows a phylogenetic signal, which rejected the lambda 0 tree.
#if at last step p < 0.05 -- Rejects no signal model, trait has phylogenetic signal
#--------------------------------------------------

library(geiger)

tree=tree2

#we rescale the tree with lambda=0 to obtain a starlike tree

tree0 <- rescale(tree, model = "lambda", 0) 


# load trait data
traits=CheckM3_annotated[,'SS']
names(traits)=CheckM3_annotated[,'Bin_Id']

# lambda 1 and lambda 0 trees combined with trait data
trait_geiger <- treedata(tree,traits)
trait_geiger_0 <- treedata(tree0,traits)

# Comparison of different models for trait evolution
trait_ER <- fitDiscrete(trait_geiger$phy, traits, type="discrete", model = "ER", niter = 1000) 
trait_SYM <- fitDiscrete(trait_geiger$phy, traits, type="discrete", model = "SYM", niter = 1000) 
trait_ARD <- fitDiscrete(trait_geiger$phy, traits, type="discrete", model = "ARD", niter = 1000) 

#likelihood ratio test  ER vs SYM
# p > 0.05 -- Not rejects ER-model
trait_d_ER_vs_SYM <- abs(2*(trait_ER$opt$lnL-trait_SYM$opt$lnL)) # 1.543
trait_p_value_ER_vs_SYM <- pchisq(trait_d_ER_vs_SYM, 3-1, lower.tail=FALSE) 

#likelihood ratio test  ER vs ARD 
# p > 0.05 -- Not rejects ER-model
trait_d_ER_vs_ARD <- abs(2*(trait_ER$opt$lnL-trait_ARD$opt$lnL)) # 8.855
trait_p_value_ER_vs_ARD <- pchisq(trait_d_ER_vs_ARD, 6-1, lower.tail=FALSE) # 0.114 

#--------------------------------------------------
#			MAKE YOUR DESCISION !!!
#--------------------------------------------------

#--------------------------------------------------

# IF YOU ACCEPT the ER model 
# phylogenetic signal test for trait trait, using Pagel's lambda
# p < 0.05 -- Rejects no signal model, trait has phylogenetic signal
trait_ER_lambda <- fitDiscrete(trait_geiger$phy,traits, type="discrete", model = "ER", transform = "lambda", niter = 1000)
trait_ER_lambda_0 <- fitDiscrete(trait_geiger_0$phy, traits, type="discrete", model = "ER", transform = "lambda", niter = 1000)
trait_d_our_tree_vs_lambda_0 <- abs(2*(trait_ER_lambda_0$opt$lnL-trait_ER_lambda$opt$lnL)) 
trait_p_our_tree_vs_lambda_0 <- pchisq(trait_d_our_tree_vs_lambda_0, 1, lower.tail=FALSE) 

#--------------------------------------------------

# IF YOU ACCEPT the SYM model 
# phylogenetic signal test for trait trait, using Pagel's lambda
# p < 0.05 -- Rejects no signal model, trait has phylogenetic signal
trait_SYM_lambda <- fitDiscrete(trait_geiger$phy,traits, type="discrete", model = "SYM", transform = "lambda", niter = 1000)
trait_SYM_lambda_0 <- fitDiscrete(trait_geiger_0$phy, traits, type="discrete", model = "SYM", transform = "lambda", niter = 1000)
trait_d_our_tree_vs_lambda_0 <- abs(2*(trait_SYM_lambda_0$opt$lnL-trait_SYM_lambda$opt$lnL)) 
trait_p_our_tree_vs_lambda_0 <- pchisq(trait_d_our_tree_vs_lambda_0, 1, lower.tail=FALSE) 






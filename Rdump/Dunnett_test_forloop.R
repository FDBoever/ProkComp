#### Dunnett's test of all-against-one comparison. This is preferred post-hoc test when comparing all treatments
#### versus a single control.

#### I've used https://stats.stackexchange.com/questions/83116/dunnetts-test-in-r-returning-different-values-each-time
#### for the multcomp aov and Dunnett's test routine
####

##### NOTE. Theuse of sink() below doesn't work when file is sourced in RStudio using, source("Dunnett_test.R")
##### Need to source file and include 'print.eval = T', i.e.  source("Dunnett_test.R", print.eval = T)

library(reshape2) # use for melt
library(multcomp) # use for aov and glht

# this data set was for FAME analysis between Axenic and all other treatments for each FAME category
fa <- read.table("FAcomp.txt", header = T, sep = "\t")
fa.m <- melt(fa)

# rename colum header 
colnames(fa.m)[colnames(fa.m) == "variable"] <- "treatment"


########################################################################
#                        FOR LOOPS IN R
########################################################################


#for clarity, you can first make a vector of all unique FA classes in your dataframe 
unique.fa = unique(fa.m$FA)

# as a silly example of for loops you can run this, it merely prints the elements of the vector we created
#for (fa.class in unique.fa){
#  print(fa.class)
#}

# the following should work?!
# same rationale, the for loop runs over all elements in the unique.fa vector and computes the steps for each

sink("Dunnett_test_output.txt", append=TRUE, split = TRUE)
set.seed(20140123)	


for (fa.class in unique.fa){
  fa.sub <- subset(fa.m, FA == fa.class)
  fa.aov <- aov(value ~ treatment, fa.sub)
  summary(fa.aov)
  summary(glht(fa.aov, linfct = mcp(treatment = "Dunnett")))
  
  #in addition you could add this code, which will create an output file, and each time append the summary statistics to it
  #out1 <- capture.output(summary(fa.aov))
  #out2 <- capture.output(summary(summary(glht(fa.aov, linfct = mcp(treatment = "Dunnett")))))
  #cat(paste("\n----------\n",fa.class), "\n-- Analysis of Variance", out1,"\n-- Post hoc test",out2, file="SummaryStatistics.txt", sep="\n", append=TRUE)
}


########################
########################



# Subsetting data based on FA category
fa14 <- subset(fa.m, FA == "14:00")
fa16 <- subset(fa.m, FA == "16:00")
fa18 <- subset(fa.m, FA == "18:00")
fa16.1.7 <- subset(fa.m, FA == "16:1(n-7)")
fa16.1.9 <- subset(fa.m, FA == "16:1(n-9)")
fa18.1.7 <- subset(fa.m, FA == "18:1(n-7)")
fa18.1.9 <- subset(fa.m, FA == "18:1(n-9)")
fa18.2.6 <- subset(fa.m, FA == "18:2(n-6)")
fa18.3.6 <- subset(fa.m, FA == "18:3(n-6)")
fa20.3.6 <- subset(fa.m, FA == "20:3(n-6)")
fa20.4.6 <- subset(fa.m, FA == "20:4(n-6)")
fa20.5.3 <- subset(fa.m, FA == "20:5(n-3)")

# Analysis of variance using multcomp
set.seed(20140123)	
aov14 <- aov(value ~ treatment, fa14)
aov16 <- aov(value ~ treatment, fa16)
aov18 <- aov(value ~ treatment, fa18)
aov16.1.7 <- aov(value ~ treatment, fa16.1.7)
aov16.1.9 <- aov(value ~ treatment, fa16.1.9)
aov18.1.7 <- aov(value ~ treatment, fa18.1.7)
aov18.1.9 <- aov(value ~ treatment, fa18.1.9)
aov18.2.6 <- aov(value ~ treatment, fa18.2.6)
aov18.3.6 <- aov(value ~ treatment, fa18.3.6)
aov20.3.6 <- aov(value ~ treatment, fa20.3.6)
aov20.4.6 <- aov(value ~ treatment, fa20.4.6)
aov20.5.3 <- aov(value ~ treatment, fa20.5.3)

# now do Dunnett's test and use  to write following output to file
# The following use of sink() doesn't work when file is sourced in RStudio e.g. source("Dunnett_test.R")
# Need to source file like this, source("Dunnett_test.R", print.eval = T)
sink("Dunnett_test_output.txt", append=TRUE, split = TRUE)
summary(aov14)
summary(glht(aov14, linfct = mcp(treatment = "Dunnett")))
summary(aov16)
summary(glht(aov16, linfct = mcp(treatment = "Dunnett")))
summary(aov18)
summary(glht(aov18, linfct = mcp(treatment = "Dunnett")))
summary(aov16.1.7)
summary(glht(aov16.1.7, linfct = mcp(treatment = "Dunnett")))
summary(aov16.1.9)
summary(glht(aov16.1.9, linfct = mcp(treatment = "Dunnett")))
summary(aov18.1.7)
summary(glht(aov18.1.7, linfct = mcp(treatment = "Dunnett")))
summary(aov18.1.9)
summary(glht(aov18.1.9, linfct = mcp(treatment = "Dunnett")))
summary(aov18.2.6)
summary(glht(aov18.2.6, linfct = mcp(treatment = "Dunnett")))
summary(aov18.3.6)
summary(glht(aov18.3.6, linfct = mcp(treatment = "Dunnett")))
summary(aov20.3.6)
summary(glht(aov20.3.6, linfct = mcp(treatment = "Dunnett")))
summary(aov20.4.6)
summary(glht(aov20.4.6, linfct = mcp(treatment = "Dunnett")))
summary(aov20.5.3)
summary(glht(aov20.5.3, linfct = mcp(treatment = "Dunnett")))
sink()



# This script is designed to visualise domain architecture of user specified gene(s) and domains.
# It will annotate the domains as well as depict the region, with user specificied data
# for now it plots everything on different lines, which allows selection and inspection
# at the current state, you may want to post-process the figure in Illustrator/InkScape, to overlay the domains with the genes of interest, and color at will. 
 
# This can serve as an example of how to build the dataframe
# this information can be obtained via, HHpred, PfamScan, Phyre2, ExPASy-PROSITE or other domain finders

###########################
#library(devtools)
#install_github("trinker/plotflow")

library(gdata)
library(ggplot2)
library(plotflow)
library(dplyr)
library(RColorBrewer)
###########################

 
 #CAUTION!
 # MAKE SURE YOU HAVE THE GENES OF INTEREST AS RECORDS IN THE DATAFRAME
 # Start=0 , end =nr of AA in total sequence
 
 
dfDomains = data.frame('Domain'=c('UPF0257','1gene','KOG4659','6FB3_B','4O9X_A','5KIS_B','COG3209','PF14436.5','gene2','LysM','5JCE_B'),'Type'=c('UPF0257','1gene','KOG4659','6FB3_B','4O9X_A','5KIS_B','COG3209','PF14436.5','gene2','LysM','5JCE_B'),'Start'=c(49,0,567,486,698,799,569, 1524,0,3, 82),'Stop'=c(211 ,1616 ,1460,1456,1475,1462,1450, 1615,702,48, 139))


plot =ggplot(data=dfDomains, aes(x = Type, ymin = Start, ymax = Stop, colour = Domain))
plot + 
  geom_linerange(size = 10) + 
  coord_flip() +
  geom_text(aes(y = Start, x = Type, label = Domain), hjust = 0, vjust = -2, size = 4) + 
  geom_text(aes(y = Start, x = Type, label = Start), hjust = 0, vjust = 3, size = 3) + 
  geom_text(aes(y = Stop, x = Type, label = Stop), hjust = 1, vjust = 3, size = 3) +
  theme_bw() +
  labs(y = NULL, x = NULL, title = "FDB33_02831") +
  theme(axis.line = element_blank(), 
  panel.grid.major =  element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background =  element_blank(), 
  axis.ticks = element_blank(), 
  legend.position ="none", 
  axis.text = element_blank())+scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(dim(dfDomains)[1]))
  

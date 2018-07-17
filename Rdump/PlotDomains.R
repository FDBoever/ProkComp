library(gdata)
library(ggplot2)
library(plotflow)
library(dplyr)

#library(devtools)
#install_github("trinker/plotflow")


dfDomains = data.frame('Domain'=c('UPF0257','gene','KOG4659','6FB3_B','4O9X_A','5KIS_B','COG3209','PF14436.5','gene2','LysM','5JCE_B'),'Type'=c('UPF0257','gene','KOG4659','6FB3_B','4O9X_A','5KIS_B','COG3209','PF14436.5','gene2','LysM','5JCE_B'),'Start'=c(49,0,567,486,698,799,569, 1524,0,3, 82),'Stop'=c(211 ,1616 ,1460,1456,1475,1462,1450, 1615,702,48, 139))


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
  axis.text = element_blank())+ylim(c(0,1616))


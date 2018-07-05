library(ggtree)
library(ape)
library(data.table)
library(RColorBrewer)
library(phytools)
library(ggplot2)
###########################
alkB_cytoscape = read.csv('~/Downloads/alkB_full_combined.csv',stringsAsFactors=FALSE)


alkB_tree = read.tree('~/Genomics/alkB_fasttree.tre')

alkB_tree=midpoint.root(alkB_tree)
#atpD_tree =midpoint.root(atpD_tree)

alkB_tree$tip.label = sub("\\_.*", "", alkB_tree$tip.label)
alkB_tree$tip.label = sub("tr|", "", alkB_tree$tip.label)

#alkB_tree$tip.label  = as.character(unique(unlist(strsplit(alkB_tree$tip.label, "[| ]")))

correctedTipLabels = strsplit(alkB_tree$tip.label, "[| ]")
correctedTipLabels = unlist(lapply(correctedTipLabels, tail, 1))
alkB_tree$tip.label = correctedTipLabels

rownames(alkB_cytoscape) = alkB_cytoscape$name
ralkB = alkB_cytoscape
ralkB$name = as.character(ralkB$name)
#ralkB = ralkB[as.character(alkB_tree$tip.label), ]
#ralkB = ralkB[grepl(paste(alkB_tree$tip.label,collapse="|"),ralkB$name), ]


ralkB = alkB_cytoscape[grepl(paste(unique(unlist(strsplit(alkB_tree$tip.label, "[| ]"))),collapse="|"), alkB_cytoscape $name), ]
rownames(ralkB) = ralkB$name
ralkB = ralkB[alkB_tree $tip.label,]


#write the following, could serve as supplementary information
alkB_cytoscape[alkB_cytoscape $PFAM == 'PF00301|PF00487',]
alkB_cytoscape[alkB_cytoscape $Genus == 'Marinobacter',]


###########
ralkB$Genera <- ralkB $Genus

ralkB$Genera[ralkB$Genus!="Marinobacter"] <- NA
ralkB$Genera[ralkB$Genus=="Alcanivorax"] <- "Alcanicorax"
ralkB$Genera[ralkB$Genus=="Pseudomonas"] <- "Pseudomonas"
ralkB$Genera[ralkB$Genus=="Acinetobacter"] <- "Actinetobacter"
ralkB$Genera[ralkB$Genus=="Rhodococcus"] <- "Rhodococcus"
ralkB$Genera[ralkB$Genus=="Mycobacterium"] <- "Mycobacterium"
ralkB$Genera[ralkB$Genus=="Burkholderia"] <- "Burkholderia"
ralkB$Genera[ralkB$Genus=="Micromonospora"] <- "Micromonospora"
ralkB$Genera[ralkB$Genus=="Psychrobacter"] <- "Psychrobacter"
ralkB$Genera[ralkB$Genus=="Streptomyces"] <- "Streptomyces"





tree = alkB_tree

my_info <- data.table(tip_lbs = ralkB$name,
                      groupA = ralkB$Genera,
                      groupB = ralkB$PFAM,
                      val = ralkB$Sequence.Length)

grA <- split(my_info$tip_lbs, my_info$groupA)
grA

tree_grA <- ggtree::groupOTU(tree, grA)
str(tree_grA)


levels(attributes(tree_grA)$group)[1] <- levels(attributes(tree_grA)$group)[2]

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#fcolors = c("#999999","#629363","#3A85A8","#C4625D","#E41A1C","#FFC81D","#BF862B","#EB7AA9")
#fcolors = getPalette(unique(my_info$groupB))
fcolors = getPalette(20)
#A_fcolors = getPalette(unique(my_info$groupA))
fcolors=c('#3D505E','#76E6CA','#E59535','#DB0000','#E6999A','#74C94C','#969696','gray45','gray85','#FFFD5D','#BFFFFE','#AD23C4','#796A33')

tree_plot <- 
  ggtree(tr = tree_grA, 
         # color by group attribute, check str(tree_grA)
         mapping = aes(color = group), 
         layout  = 'circular',
         # set line thikness
         size = 0.4) + scale_color_manual(name = 'Group A',
                     values = fcolors)


tree_plot <- tree_plot %<+% my_info 

tree_dt <- data.table(tree_plot$data)
head(tree_dt)

# select only the tip labels and order by coord y
tree_dt <- tree_dt[isTip == TRUE][order(y)]

# Make table with y cords for each groupB in consecutive order; 
# this helps for drawing & labeling segments.
# Note the usage of "rleid" function, which is a grouping ID generator,
# needed because consecutive rows of an identical reoccurring group must form a unique group.
coord_groups <- tree_dt[, .(y1 = y[1],
                            y2 = y[.N],
                            angle = mean(angle),
                            n = .N), # optional - helps with counting
                        by = .(groupB, 
                               id_gr = rleid(groupB, 
                                             prefix = "grp"))]
coord_groups

# Compute the middle y - will be used for placing the group label;
# similarly the mean angle was computed above already. 
coord_groups[, y_mid := rowMeans(.SD), .SDcols = c("y1", "y2")]

# For one-record groups where y1=y2, adjust their y coordinates so that a segment gets drawn;
# If not, no segment get drawn for such cases.
# To force a segment add and subtract a small amount (try & error until seems ok). 
# Prefer a smaller values since with bigger ones you risk to exaggerate the segments.
coord_groups[, y1_adj := ifelse(y1 == y2, y1 - 0.1, y1)]
coord_groups[, y2_adj := ifelse(y1 == y2, y2 + 0.1, y2)]

# Labels need angle adjustment for cases between 90 and 270 dg
coord_groups[, angle_adj := ifelse(angle %between% c(90, 180), 
                                   yes = angle + 180,
                                   no = ifelse(angle > 180 & angle <= 270,
                                               yes = angle - 180,
                                               no = angle))]

# Labels with angles between 90 and 270 dg
# need change of horizontal adjustment argument from 0 to 1.
coord_groups[, hjust_adj := ifelse(angle %between% c(90, 270), yes = 1L, no = 0L)]

# if needed, coloring could be binary 
# coord_groups[, col := ifelse(.I%%2, 0.5, 1)]
coord_groups
colourCount = length(unique(coord_groups$groupB))


# Define variable to control x coordinate of segments & labels
my_x <- max(tree_dt$x) + 0.05

tree_labeled <- 
  tree_plot + 
  	geom_segment(data = coord_groups,aes(x = my_x, y = y1_adj, xend = my_x, yend = y2_adj,color = groupB),lineend = "butt",size = 3) +
 	geom_text(data = coord_groups,
aes(x = my_x,y = y_mid,angle = angle_adj,hjust = hjust_adj,label = groupB),vjust = 0.5, size  = 2.5,nudge_x = 0.05, color = "black") +
  theme(
    text = element_text(size = 10, family = "sans"),
    legend.justification = c(1,0),
    legend.position = c(1.1, 0.05),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.key.height = unit(4, "mm"),
    plot.margin = unit(c(t = 0.2, r = 1.3, b = -0.35, l = -0.2), "cm")
  )
tree_labeled





coord_groups2 <- tree_dt[, .(y1 = y[1],y2 = y[.N],angle = mean(angle),n = .N), # optional - helps with counting
                        by = .(groupA,id_gr = rleid(groupA,prefix = "grp"))]

coord_groups2[, y_mid := rowMeans(.SD), .SDcols = c("y1", "y2")]
coord_groups2[, y1_adj := ifelse(y1 == y2, y1 - 0.1, y1)]
coord_groups2[, y2_adj := ifelse(y1 == y2, y2 + 0.1, y2)]
coord_groups2[, angle_adj := ifelse(angle %between% c(90, 180), yes = angle + 180,no = ifelse(angle > 180 & angle <= 270,yes = angle - 180,no = angle))]
coord_groups2[, hjust_adj := ifelse(angle %between% c(90, 270), yes = 1L, no = 0L)]

colourCountA = length(unique(coord_groups$groupA))






# Define variable to control the x coordinates of bars (segments)
my_factor <- 0.0003
x_base <- max(tree_dt$x) + abs(min(tree_dt$val, na.rm = TRUE))*my_factor + 0.02

# Define variable to control the x coordinates of segments & labels
my_x <- x_base + max(tree_dt$val, na.rm = TRUE)*my_factor + 0.05

# Need to add a value (usually a small amount) to `ymax` 
# to force a complete circle (otherwise a stripe of white remains).
# This value needs to be just right:
# -  not too big because affects the visual which results in strange angles on labels;
# -  not too small because otherwise the remaining strip of white is not completely filled.
# Is a matter of try and error until you get it right. 
# Also, better to be biased towards smaller values since big value can lead to 
# big displacements when drawing the tree.
fill_in_value <- 0.2

# Categorize values - binary
tree_dt[, categ := ifelse(val <= 0, "neg", "pos")]

# Plot the tree with circular barplot
tree_bars <- 
  tree_plot +
  # Add a disc to plot bars on top of it
  geom_rect(data = tree_dt,
            aes(xmin = x_base + min(val, na.rm = TRUE)*my_factor,
                ymin = 0,
                xmax = x_base + max(val, na.rm = TRUE)*my_factor,
                # Add also fill_in_value to `ymax` to force a complete circle
                ymax = max(y) + fill_in_value), 
            color = NA, # set NA so to avoid coloring borders
            fill = "#deebf7", # or try "#8da0cb"
            alpha = 0.1) +
  # Add bars for each tip label
  geom_rect(data = tree_dt,
            aes(xmin = x_base,
                ymin = y - 0.5,
                xmax = x_base + val*my_factor,
                ymax = y + 0.5,
                fill = categ),
            # no borders
            color = NA) +
  # Fill the bars
  scale_fill_manual(name   = 'Some variable',
                    breaks = c("neg", "pos"),
                    values = c("neg" = "#fc8d62",
                               "pos" = "darkseagreen4"),
                    labels = c("value1", "value2")) +
  # Add a ring (serves as reference for bars)
  geom_segment(data = tree_dt,
               aes(x = x_base, 
                   y = 0, 
                   xend = x_base, 
                   # Add also fill_in_value to `yend` to force a complete circle
                   yend = max(y) + fill_in_value),
               color = "gray50",
               # linetype = "dashed", 
               size = 0.25) +
  # Add line segments for each group.
  geom_segment(data = coord_groups,
               aes(x = my_x, 
                   y = y1_adj, 
                   xend = my_x, 
                   yend = y2_adj, color = groupB),
              
               lineend = "butt",
               size = 3) +
               	geom_segment(data = coord_groups2,aes(x = my_x+0.07, y = y1_adj, xend = my_x+ +0.07, yend = y2_adj,color = groupA),lineend = "butt",size = 3)+
  # Add text group labels at the middle of each segment.
 geom_text(data = coord_groups2,
aes(x = my_x,y = y_mid,angle = angle_adj,hjust = hjust_adj,label = groupA),vjust = 0.5, size  = 2.5,nudge_x = 0.15, color = "black") +
  theme(
    # Set font size & family - affects legend only 
    # "sans" = "Arial" and is the default on Windows OS; check windowsFonts()
    text = element_text(size = 10, family = "sans"),
    # Grab bottom-right (x=1, y=0) legend corner 
    legend.justification = c(1,0),
    # and position it in the bottom-right plot area.
    legend.position = c(1.25, 0.1),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    # Set height of legend items (keys). This affects spacing between them as well.
    legend.key.height = unit(4, "mm"),
    # Set margin around entire plot.
    plot.margin = unit(c(t = 0.4, r = 2.6, b = -0.6, l = 0), "cm")
  )

tree_bars











tree_labeled2 <- 
  tree_plot + 
  	geom_segment(data = coord_groups,aes(x = my_x, y = y1_adj, xend = my_x, yend = y2_adj,color = groupB),lineend = "butt",size = 3) +
	geom_segment(data = coord_groups2,aes(x = my_x+0.07, y = y1_adj, xend = my_x+ +0.07, yend = y2_adj,color = groupA),lineend = "butt",size = 3) +

 	geom_text(data = coord_groups2,
aes(x = my_x,y = y_mid,angle = angle_adj,hjust = hjust_adj,label = groupA),vjust = 0.5, size  = 2.5,nudge_x = 0.15, color = "black") +
  theme(
    text = element_text(size = 10, family = "sans"),
    legend.justification = c(1,0),
    legend.position = c(1.1, 0.05),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.key.height = unit(4, "mm"),
    plot.margin = unit(c(t = 0.2, r = 1.3, b = -0.35, l = -0.2), "cm")
  )
tree_labeled2

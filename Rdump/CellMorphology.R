library("EBImage")

#f = system.file("/Volumes/NO NAME/august_018/SNAP-141208-0092.jpg", package="EBImage")
img = readImage("/Volumes/NO NAME/august_018/SNAP-141208-0092.jpg")
img = readImage("/Volumes/NO NAME/august_018/SNAP-135707-0065.jpg")
#img = readImage('/Volumes/NO NAME/august_018/SNAP-122134-0046.jpg')

display(img, method="raster")
/Volumes/NO NAME/august_018/SNAP-135707-0065.jpg
hist(img)
range(img)


img.gray <-channel(img,"gray")
display(img.gray, method="raster")

threshold = otsu(img.gray)
threshold

nuc_th = combine( mapply(function(frame, th) frame > th, getFrames(img.gray), threshold, SIMPLIFY=FALSE) )


disc = makeBrush(31, "disc")
disc = disc / sum(disc)
offset = 0.05
nuc_bg = filter2( img.gray, disc )
nuc_th = img.gray > nuc_bg + offset


display(nuc_th, all=TRUE,method="raster")





nmask = watershed( distmap(nuc_th), 2 )
display(colorLabels(nmask), all=TRUE,method="raster")



nmask = thresh(img.gray, w=10, h=10, offset=0.05)
nmask = opening(nmask, makeBrush(5, shape='disc'))
nmask = fillHull(nmask)
nmask = bwlabel(nmask)

#nmask=max(nmask)-nmask

display(nmask, all=TRUE,method="raster")
display(max(nuc_th)-nuc_th, all=TRUE,method="raster")

#st = stackObjects(nmask, img)
st = stackObjects(nmask, nmask)
st2 = stackObjects(nmask, img)

#st = stackObjects(max(nuc_th)-nuc_th, max(nuc_th)-nuc_th)


display(st, all = TRUE,method="raster")
display(st2, all = TRUE,method="raster")

st_sel=st@.Data[,,as.numeric(rownames(df.img[df.img$s.area>300,]))]
display(st_sel, all = TRUE,method="raster")


writeImage(st, "Alteromonas_test.jpeg", quality = 100)


df.img = data.frame(computeFeatures.shape(nmask))
testhist2 = hist(df.img$s.area)
plot(df.img$s.area,df.img$s.perimeter)
plot(df.img)





#

install.packages("Momocs", dependencies=TRUE, repos='http://cran.rstudio.com/')
library(Momocs)

library(magick)



lf <- list.files('~/momocs_test', full.names=TRUE)

#make images negative

for(img_f in lf){
	imgf = readImage(img_f)
	image_write(image_negate(imgf), path = img_f, format = "jpg")

}


coo <- import_jpg(lf[1:100])
out.coo = Out(coo)





stack(out.coo, title="Raw Outlines")

# 1. smoothing -> noise reduction from digitalization process #
shape.coo.adj1<-coo_smooth(out.coo, 5)
stack(shape.coo.adj1, title="Smoothed")

# 2. centering -> all centroid moved to same coordinate
shape.coo.adj2<-coo_center(shape.coo.adj1)
stack(shape.coo.adj2, title ="Centred")

# 3. re-scaling -> size-correction
shape.coo.adj3<-coo_scale(shape.coo.adj2)
stack(shape.coo.adj3, title ="Scaled")

shape.coo.adj4<-coo_slidedirection(shape.coo.adj3,"S")
stack(shape.coo.adj4, title="Normalised")

# 4. Sampling Pseudo-Landmarks "Points sampled"
# Initially 1000: but now less coordinates than n pseudo-landmarks due to smaller .jpg size
# reduces file size from 6 to 3.4Mb
shape.coo.adj5<-coo_sample(shape.coo.adj4, 40)
stack(shape.coo.adj5, title="Sample Points", pch=20, cex=1.3, points = TRUE)


shape.coo.adj6<-fgProcrustes(shape.coo.adj5, tol=0.00000000000000000000000000001)
stack(shape.coo.adj6, title="Aligned")

shape.coo.adj7<-fgProcrustes(shape.coo.adj6)
stack(shape.coo.adj7, title="Aligned")

fgsProcrustes(Ldk(shape.coo.adj5))
pProcrustes(shape.coo.adj5, coo[4])


plot(coo[1])

outline <- Out(raw_data$outline)
ef <- efourier(out.coo)

landmarks <- Ldk(out.coo)

fg <- fgProcrustes(Ldk(shape.coo.adj5))

stack(fg, title="Aligned")







mussel.group <- lf_structure(lf, names = c("group","ID"), split = "_", trim.extension = TRUE)

# Coo-Object = Collection of coordinates; 
# along with grouping factors e.g. "ID", "Sampling", "Depth" #
mussel.coo<-Out(mussel.out,fac=mussel.group[1:3])
mussel.coo



### -----      II. OUTLINE INSPECTION    ----- ####
### > Panel Plot ----
# all info
panel(mussel.coo, names=TRUE)
panel(LS.coo, names=TRUE)

# as factor "Depth" (names="Sampling" possible)
panel(mussel.coo, fac="Depth", palette = col_heat, names="ID", main="Mussel Shell Shapes per Depth")
panel(LS.coo, fac="Depth", palette = col_heat, names="ID", main="Shell Shapes per Depth - Loch Spelve")


### Center Image(s) - "Try-out"
# By way of example using a single shape:
# test.center<-mussel.coo[4]
# par(mfrow=c(1,2))
# coo_plot(test.center, main = "004_0_1")
# coo_plot(coo_center(test.center),main="centered 004_0_1")

## or, using elliptic fourier transforms, fitting x and y coordinates separately #
# coo_oscillo(mussel.coo[4],"efourier")

    
###      III. OUTLINE ADJUSTMENT & NORMALIZATION    ####
# i.e. for following shape coordinates to be invariant to outline size, 
# rotation and position
# each step can be plotted with stack(mussel.coo.adj, title="...", borders="put a colour here")
stack(test, title="Raw Outlines (n = 431)")

# 1. smoothing -> noise reduction from digitalization process #
mussel.coo.adj1<-coo_smooth(test, 5)
stack(mussel.coo.adj1, title="Smoothed")

# 2. centering -> all centroid moved to same coordinate
mussel.coo.adj2<-coo_center(mussel.coo.adj1)
stack(mussel.coo.adj2, title ="Centred")

# 3. re-scaling -> size-correction
mussel.coo.adj3<-coo_scale(mussel.coo.adj2)
stack(mussel.coo.adj3, title ="Scaled")

mussel.coo.adj4<-coo_slidedirection(mussel.coo.adj3,"S")
stack(mussel.coo.adj4, title="Normalised")

# 4. Sampling Pseudo-Landmarks "Points sampled"
# Initially 1000: but now less coordinates than n pseudo-landmarks due to smaller .jpg size
# reduces file size from 6 to 3.4Mb
mussel.coo.adj5<-coo_sample(mussel.coo.adj4, 500)
stack(mussel.coo.adj5, title="Sample Points", pch=20, cex=1.3, points = TRUE)

LS.coo.adj<-coo_sample(LS.coo.adj, 500)
stack(LS.coo.adj, title="LS - Sample Points", pch=20, cex=1.3, points = TRUE)

mussel.coo.adj6<-fgProcrustes(mussel.coo.adj4)
stack(mussel.coo.adj6, title="Aligned")

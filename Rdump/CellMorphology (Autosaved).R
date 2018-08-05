library("EBImage")

#f = system.file("/Volumes/NO NAME/august_018/SNAP-141208-0092.jpg", package="EBImage")
img = readImage("/Volumes/NO NAME/august_018/SNAP-141208-0092.jpg")
img = readImage("/Volumes/NO NAME/august_018/SNAP-135707-0065.jpg")
#img = readImage('/Volumes/NO NAME/august_018/SNAP-122134-0046.jpg')

#img = readImage('~/Downloads/s/1.15.jpg')

~/Downloads/s/1.15.jpg
display(img, method="raster")
/Volumes/NO NAME/august_018/SNAP-135707-0065.jpg
hist(img)
range(img)

#img=max(img)-img
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


shape.coo.xaln <- coo_alignxax(shape.coo.adj1)
stack(shape.coo.xaln, title="Sample Points", pch=20, cex=1.3, points = TRUE)




# 2. centering -> all centroid moved to same coordinate
#shape.coo.adj2<-coo_center(shape.coo.adj1)

shape.coo.adj2<-coo_center(shape.coo.xaln)
stack(shape.coo.adj2, title ="Centred")

# 3. re-scaling -> size-correction
#shape.coo.adj3<-coo_scale(shape.coo.adj2)
#stack(shape.coo.adj3, title ="Scaled")
shape.coo.adj3 = shape.coo.adj2
#shape.coo.adj4<-coo_slidedirection(shape.coo.adj3,"S")

shape.coo.adj4<-coo_slidedirection(shape.coo.adj3)
stack(shape.coo.adj4, title="Normalised")

# 4. Sampling Pseudo-Landmarks "Points sampled"
# Initially 1000: but now less coordinates than n pseudo-landmarks due to smaller .jpg size
# reduces file size from 6 to 3.4Mb



shape.coo.adj5<-coo_sample(shape.coo.adj4, 20)
stack(shape.coo.adj5, title="Sample Points", pch=20, cex=1.3, points = TRUE)



shape.coo.adj6<-fgProcrustes(shape.coo.xaln, tol=0.00000000000000000000000000001)

shape.coo.adj7 = coo_alignxax(shape.coo.adj6)
shape.coo.adj8 = coo_slidedirection(shape.coo.adj7)
stack(shape.coo.adj8)





stack(coo_alignxax(shape.coo.adj6))

stack(shape.coo.adj6, title="Aligned")


### SAA 2015 -- Open Methods Section (April 2015) ###
### Andrew Bevan (University College London, a.bevan_at_ucl.ac.uk) ###

## See also Bevan, A. 2012. Spatial methods for analysing large-scale artefact inventories, Antiquity 86.332: 492-506.##

##########################################

## 1. Set-Up ##

#Set a working directory
setwd("~/Desktop/Bevan") #MacOSX

# Add libraries
library(rgdal) # spatial data interoperability
library(maptools) # various methods for spatial data
library(spatstat) # point process models and other utilities

# Add a couple of custom functions for scalebar, north arrow etc.
source("utilities.R")

##########################################

## 2. Basic Data Manipulation ##

# Load basic geogrpahy for UK West Country (from http://gadm.org, generalised)
wc <- readOGR("data/shp/wc", "wc")

# Read in Iron Age coin data (from the Celtic Coin Index, within the Portable Antiquities Scheme; http://finds.org.uk) and convert to spatial data
iacoins <- read.csv("data/csv/iacoins.csv", header=TRUE)
coordinates(iacoins) <- ~Easting + Northing
proj4string(iacoins) <- CRS(proj4string(wc))

# Example of dataset re-location to obfuscate real locations without loosing precision.
set.seed(101)
myshift <- c(10000*runif(1),10000*runif(1))
wc <- elide(wc, shift=myshift)
iacoins <- elide(iacoins, shift=myshift)

# Create a convenient bounding box polygon
wcbox <- bbox(wc)
wcbox <- Polygon(cbind(c(wcbox[1,1], wcbox[1,1], wcbox[1,2], wcbox[1,2], wcbox[1,1]), c(wcbox[2,1], wcbox[2,2], wcbox[2,2], wcbox[2,1], wcbox[2,1])))
wcbox <- Polygons(list(wcbox), "1")
wcbox <- SpatialPolygons(list(wcbox), 1:1, proj4string=CRS(proj4string(wc)))

# Jitter coin findspot coordinates by a metre to remove superimposed finds
set.seed(101)
xys <- coordinates(iacoins)
xys[ ,1] <-  xys[ ,1] + runif(xys[ ,1], -1,1)
xys[ ,2] <-  xys[ ,2] + runif(xys[ ,2], -1,1)
iacoins <- iacoins@data
iacoins$Easting <- xys[,1]
iacoins$Northing <- xys[,2]
coordinates(iacoins) <- ~Easting + Northing
proj4string(iacoins) <- CRS(proj4string(iacoins))

##########################################

## 3. Regional Analysis ##

# Clean up metal attributes
golddenoms <- c("Stater (gold)","Quarter stater (gold)")
iacoins$Mat[iacoins$PrimaryMaterial=="Gold" | iacoins$Denomination %in% golddenoms] <- "Au"
iacoins$Mat[is.na(iacoins$Mat)] <- "Other"

# Extract those coins associated with the so-called Dobunni tribal area, and all others.
dobunni <- iacoins[iacoins$Tribe=="Dobunni",]
otherIA <- iacoins[iacoins$Tribe!="Dobunni",]

# Euclidean mean centre of Dobunni-type coins
dmc <- SpatialPoints(data.frame(X=mean(coordinates(dobunni)[ ,1]), Y=mean(coordinates(dobunni)[ ,2])), proj4string=CRS(proj4string(dobunni)))

# Location and size of north arrow, scalebar
xpos <- 250000
ypos <- 110000
scalesize <- 50000

# Plot
plot(wc, col="grey75", border="grey75")
plot(wcbox, add=TRUE)
plot(otherIA, pch=19, cex=0.2, col="grey25", add=TRUE)
plot(dobunni, pch=19, cex=0.5, col="red", add=TRUE)
plot(dmc, pch=15, cex=1.5, col="yellow", add=TRUE)
legend("topleft",legend=c("Dobunni-type","Other coins","Mean centre"), col=c("red","grey25","yellow"), pch=c(19,19,15), pt.cex=c(0.5,0.2,1.5), bty='n', cex=0.8, inset=c(0.12,0.05))
northArrow(xpos, ypos, scalesize / 5, lwd=0.5)
scaleBar(xpos + (0.5 * scalesize), ypos - (0.5 * scalesize), scalesize, scalesize / 6, lwd=0.5, col="white")
text(xpos, ypos - (0.7*scalesize), labels="50 km", cex=0.75)

# Convert spatstat objects
dobC <- as.ppp(coordinates(dobunni),as.owin(wc))
othC <- as.ppp(coordinates(otherIA),as.owin(wc))

# Kernel density surface with a 1sd=7.5km Gaussian kernel and 500m cell size.
dobdens <- density(dobC, sigma=7500, edge=TRUE, eps=500)
dobdens[as.matrix(dobdens) < 0] <- 0

# Create a marked point pattern for gold Dobunni-type coins and those of other metal (silver or silver wash).
mDobAu <- as.ppp(coordinates(dobunni[dobunni$Mat == 'Au',]),as.owin(wc))
marks(mDobAu) <- as.factor("Au")
mDobOth <- as.ppp(coordinates(dobunni[dobunni$Mat == 'Other',]),as.owin(wc))
marks(mDobOth) <- as.factor("Other")
mDobC <- superimpose(mDobOth,mDobAu)
# Create a relative risk surface of gold versus other metal
rrDobC <- relrisk(mDobC, sigma=7500, edge=TRUE, eps=500)
# Remove mapping of this surface in areas of extremley low overall density (arbitrary eyeballed threshold).
denscutoff <- 0.000000004
rrDobC[as.matrix(dobdens)< denscutoff] <- NA

# Colour ramp
tc <- colourmap(rev(heat.colors(10)), breaks=c(seq(0,1,0.1)))

# Plot
dev.new(width=8, height=4)
par(mfrow=c(1,2))
par(mar=c(0, 0, 0, 0))
plot(wcbox)
plot(wc, col="grey65", border=NA, add=TRUE)
plot(dobunni[dobunni$Mat == 'Other',], pch=15, cex=0.5, col="cyan", add=TRUE)
plot(dobunni[dobunni$Mat == 'Au',], pch=19, cex=0.4, col="brown3", add=TRUE)
plot(dmc, pch=15, cex=0.6, add=TRUE)
plot(wcbox, add=TRUE)
northArrow(xpos, ypos, scalesize / 5, lwd=0.5)
scaleBar(xpos + (0.5 * scalesize), ypos - (0.5 * scalesize), scalesize, scalesize / 6, lwd=0.5, col="white")
text(xpos, ypos - (0.7*scalesize), labels="50 km", cex=0.75)
title("Gold (red) and other 'Dobunnni'-type coins", line=-2.3, cex.main=0.8)
plot(wcbox)
plot(wc, col="grey65", border=NA, add=TRUE)
plot(rrDobC, col=tc, add=TRUE)
plot(tc, vertical=T, add=T, xlim=c(460000,470000), ylim=c(190000,270000), cex.axis=0.7)
plot(dobunni[dobunni$Mat == 'Other',], pch=15, cex=0.5, col="cyan", add=TRUE)
plot(dobunni[dobunni$Mat == 'Au',], pch=19, cex=0.4, col="brown3", add=TRUE)
plot(dmc, pch=15, cex=0.6, add=TRUE)
plot(wcbox, add=TRUE)
title("Relative risk of gold coins", line=-2.3, cex.main=0.8)

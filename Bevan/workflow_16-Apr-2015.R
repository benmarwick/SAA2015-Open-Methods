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

# Example of dataset re-location
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

# Clean up materials
golddenoms <- c("Stater (gold)","Quarter stater (gold)")
iacoins$Mat[iacoins$PrimaryMaterial=="Gold" | iacoins$Denomination %in% golddenoms] <- "Au"
iacoins$Mat[is.na(iacoins$Mat)] <- "Other"

# Extract those coins associated with so-called Dobunni tribal area, and all others.
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

# Created a marked point pattern for gold Dobunni-type coins and those of other materials.
mDobAu <- as.ppp(coordinates(dobunni[dobunni$Mat == 'Au',]),as.owin(wc))
marks(mDobAu) <- as.factor("Au")
mDobOth <- as.ppp(coordinates(dobunni[dobunni$Mat == 'Other',]),as.owin(wc))
marks(mDobOth) <- as.factor("Other")
mDobC <- superimpose(mDobOth,mDobAu)
# Create a relative risk surface of gold versus other materials (silver or silver wash)
rrDobC <- relrisk(mDobC, sigma=7500, edge=TRUE, eps=500)
# Remove mapping of this surface in areas of extremley low overall density (arbitrary eyeballed threshold)
denscutoff <- 0.000000004
rrDobC[as.matrix(dobdens)< denscutoff] <- NA

#  Colour ramp and plot
tc <- colourmap(rev(heat.colors(10)), breaks=c(seq(0,1,0.1)))
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
plot(wcbox)
plot(wc, col="grey65", border=NA, add=TRUE)
plot(rrDobC, col=tc, add=TRUE)
plot(tc, vertical=T, add=T, xlim=c(460000,470000), ylim=c(190000,270000), cex.axis=0.7)
plot(dobunni[dobunni$Mat == 'Other',], pch=15, cex=0.5, col="cyan", add=TRUE)
plot(dobunni[dobunni$Mat == 'Au',], pch=19, cex=0.4, col="brown3", add=TRUE)
plot(dmc, pch=15, cex=0.6, add=TRUE)
plot(wcbox, add=TRUE)


# Tilman's work:
identical_windows_wrk <- function(w1,w2){
	w1x <- vertices(w1)$x
	w1y <- vertices(w1)$y
	w2x <- vertices(w2)$x
	w2y <- vertices(w2)$y
	
	if(!all(c(length(w1x)==length(w2x),length(w1y)==length(w2y)))) return(FALSE)
	
	if(any(c(sum(w1x!=w2x),sum(w1y!=w2y))>0)) return(FALSE)
	else return(TRUE)
}
	

risktol <- function(cases,controls,sigma.rrs,sigma.pooled=sigma.rrs,test="upper",zero.clip=TRUE,log.transform=TRUE,...){
	if(!identical_windows_wrk(cases$window,controls$window)) stop("cases and controls appear to have different windows")
	
	pooled.density.spill <- density(ppp(x=c(cases$x,controls$x),y=c(cases$y,controls$y),window=cases$window),sigma=sigma.pooled,...,spill=1)
	case.density.spill <- density(cases,sigma=sigma.rrs,...,spill=1)
	control.density.spill <- density(controls,sigma=sigma.rrs,...,spill=1)
	
	dummydens <- density(cases,sigma=sigma.rrs,...)
	xdatarange <- sort(rep(dummydens$xcol,length(dummydens$yrow)))
	ydatarange <- rep(dummydens$yrow,length(dummydens$xcol))
	
	case.qhz <- as.vector(case.density.spill$edg$v)
	case.qhz[!inside.owin(xdatarange,ydatarange,cases$window)] <- NA
	control.qhz <- as.vector(control.density.spill$edg$v)
	control.qhz[!inside.owin(xdatarange,ydatarange,controls$window)] <- NA
	pooled.qhz <- as.vector(pooled.density.spill$edg$v)
	pooled.qhz[!inside.owin(xdatarange,ydatarange,cases$window)] <- NA
	case.raw.v <- as.vector(case.density.spill$raw$v)
    case.raw.v[!inside.owin(xdatarange,ydatarange,cases$window)] <- NA
    case.zv <- case.raw.v/case.qhz
    
    if(zero.clip) case.zv[case.zv < 0] <- 0
    case.zv <- case.zv/sum(case.zv*dummydens$xstep*dummydens$ystep,na.rm=T)
    case.Zmat <- matrix(case.zv,length(dummydens$xcol),length(dummydens$yrow),byrow=TRUE) #?xcol v yrow
    control.raw.v <- as.vector(control.density.spill$raw$v)
    control.raw.v[!inside.owin(xdatarange,ydatarange,controls$window)] <- NA
    control.zv <- control.raw.v/control.qhz
    if(zero.clip) control.zv[control.zv < 0] <- 0
    control.zv <- control.zv/sum(control.zv*dummydens$xstep*dummydens$ystep,na.rm=T)
    control.Zmat <- matrix(control.zv,length(dummydens$xcol),length(dummydens$yrow),byrow=TRUE) #?xcol v yrow
    pooled.raw.v <- as.vector(pooled.density.spill$raw$v)
    pooled.raw.v[!inside.owin(xdatarange,ydatarange,cases$window)] <- NA
    pooled.zv <- pooled.raw.v/pooled.qhz
    if(zero.clip) pooled.zv[pooled.zv < 0] <- 0
    pooled.zv <- pooled.zv/sum(pooled.zv*dummydens$xstep*dummydens$ystep,na.rm=T)
    pooled.Zmat <- matrix(pooled.zv,length(dummydens$xcol),length(dummydens$yrow),byrow=TRUE) #?xcol v yrow
    
    rrs <- case.Zmat/control.Zmat
	if(log.transform) rrs <- log(rrs)  
  
    k2fix  <- (1/(4*pi))*as.vector(density(ppp(x=c(cases$x,controls$x),y=c(cases$y,controls$y),window=cases$window),sigma=sqrt(0.5*sigma.pooled^2),...,spill=1)$edg$v)
    RrzK <- k2fix/(as.vector(pooled.density.spill$edg$v)^2)
    denominator <- sqrt(RrzK*(cases$n^(-1)+controls$n^(-1)))/(sigma.pooled*sqrt(as.vector(t(pooled.Zmat))))
    numerator <- as.vector(t(rrs)) - !log.transform
    Zstandard <- numerator/denominator
    if(test=="upper"){
    	P <- pnorm(Zstandard,lower.tail=F)
    } else if(test=="lower"){
        P <- pnorm(Zstandard,lower.tail=T)
    } else {
        P <- 2*pnorm(abs(Zstandard),lower.tail=F)
    }
    P[!inside.owin(xdatarange,ydatarange,cases$window)] <- NA
   
	return(list(f=im(t(case.Zmat),xcol=dummydens$xcol,yrow=dummydens$yrow),
	            g=im(t(control.Zmat),xcol=dummydens$xcol,yrow=dummydens$yrow),
	            rho=im(t(rrs),xcol=dummydens$xcol,yrow=dummydens$yrow),
	            pooled=im(t(pooled.Zmat),xcol=dummydens$xcol,yrow=dummydens$yrow),
	            P=matrix(P,length(dummydens$xcol),length(dummydens$yrow),byrow=T)))
}


# plotting 2 (re-creating conditional probability, but including log density-ratio contours)
fdata <- split(mDobC)[[2]]
gdata <- split(mDobC)[[1]]
mdD.risktol <- risktol(fdata,gdata,sigma.rrs=7500,edge=TRUE,eps=500)
mdD.cp <- im(mdD.risktol$f$v/(mdD.risktol$f$v+mdD.risktol$g$v),xcol=mdD.risktol$f$xcol,yrow=mdD.risktol$f$yrow)
mdD.risktol$pooled[as.matrix(dobdens)< denscutoff] <- NA
mdD.risktol$P[as.matrix(dobdens)< denscutoff] <- NA
mdD.cp[as.matrix(dobdens)< denscutoff] <- NA

tc <- colourmap(rev(heat.colors(10)), breaks=c(seq(0,1,0.1)))
par(mar=c(0, 0, 0, 0))
plot(wcbox)
plot(wc, col="grey65", border=NA, add=TRUE)
plot(mdD.cp, col=tc, legend=FALSE, add=TRUE)
plot(dobunni[dobunni$Mat == 'Other',], pch=15, cex=0.5, col="cyan", add=TRUE)
plot(dobunni[dobunni$Mat == 'Au',], pch=19, cex=0.4, col="brown3", add=TRUE)
contour(mdD.risktol$pooled$xcol,mdD.risktol$pooled$yrow,mdD.risktol$P,levels=c(0.005,0.000001),add=T,lty=c(2,1))

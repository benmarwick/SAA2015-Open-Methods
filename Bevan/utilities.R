### SAA 2015 -- Open Methods Section (April 2015) ###
### Andrew Bevan (University College London, a.bevan_at_ucl.ac.uk) ###

##########################################

## Custom Functions ##

scaleBar <- function(x, y, width, height,...){
  #Crude scale bar for mapping
  polygon(x=c(x, (x-width), (x-width), x), y=c(y, y, (y+height), (y+height)), ...)
}

##

spEllipse <- function(xc, yc, a, b=a, angle=0, nsteps=360, proj4string=CRS("NA"), noSp=FALSE, ...){
  #Circles and ellipses, used here by northArrow()
  theta <- seq(0,(2*pi),len=nsteps)
  angle <- (90 - angle) * (pi /180)

  x <- xc + ((a * cos(theta) * cos(angle)) - (b * sin(theta) * sin(angle)))
  y <- yc + ((a * cos(theta) * sin(angle)) + (b * sin(theta) * cos(angle)))

  cpol <- cbind(x,y)
  cpol <- rbind(cpol, cpol[1,])
  if (noSp){
    return(cpol)
  } else {
    cpolp <- Polygons(list(Polygon(cpol)), ID="1")
    cpolsp <- SpatialPolygons(list(cpolp), proj4string=proj4string)
    return(cpolsp)
  }
}

##

northArrow <- function(x, y, r, type="simple", fill="black", bkgrd=NA, ...){
  #Crude north arrow for mapping
  if (type=="simple"){
    polygon(spEllipse(x,y,r, noSp=TRUE), col="white", ...)
    bl <- c((x + r * sin(210*pi/180)),(y + r * cos(210*pi/180)))
    br <- c((x + r * sin(150*pi/180)),(y + r * cos(150*pi/180)))
    polygon(x=c(x, br[1], x,  bl[1],x),y=c(y+r, br[2], (y-(0.5 * r)), bl[2],y+r), col=fill, ...)
  } else {stop("This north arrow type does not exist. Try leaving `type' blank.")}
}

##

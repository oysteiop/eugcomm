##############################################################################
#### EUGCOMM: A DATABASE OF EUGLOSSINE BEE ASSEMBLAGES ON FRAGRANCE BAITS ####
##############################################################################

rm(list=ls())

library(reshape2)
library(plyr)

# The following code demonstrates how to extract data from the R list object EUGCOMM

load("Database/R files/EUGCOMM.RData")
names(EUGCOMM)

# Extrating a single study
EUGCOMM[[14]]

allMetaData = rbind.fill(lapply(EUGCOMM, function(x) x$Metadata))

allSpeciesData = rbind.fill(lapply(EUGCOMM, function(x) x$SpeciesData))

# Replace NA with 0
for(i in 1:nrow(allSpeciesData)){
  for(j in 1:ncol(allSpeciesData)){
    if(is.na(allSpeciesData[i,j])){
      allSpeciesData[i,j]=0
    }
  }
}

# Bait data
allBaitData = rbind.fill(lapply(EUGCOMM, function(x) x$BaitData))

# Replace NA with 0
for(i in 1:nrow(allBaitData)){
  for(j in 1:ncol(allBaitData)){
    if(is.na(allBaitData[i,j])){
      allBaitData[i,j]=0
    }
  }
}

# Combine all data
AllData = data.frame(allMetaData, allBaitData, allSpeciesData)

# Species list
spList = sort(colnames(allSpeciesData[,-1]))

# Species richness
SR = rowSums(allSpeciesData[,-1]>0)
hist(SR)

# Latitudinal trend in SR
library(mgcv)
m2 = gamm(SR ~ s(lat, bs="tp"), random=list(study_area=~1), data=allMetaData)
x = seq(min(allMetaData$lat, na.rm=T), max(allMetaData$lat, na.rm=T), .1)
p = predict(m2$gam, newdata=list(lat=x), se.fit=T)

x11(height=3.5, width=8)
mat = matrix(c(1,2,3), nrow=1, byrow=T)
layout(mat=mat, widths=c(.75*.67,.25*.67,.33))
par(mar=c(4,4,1,0))
plot(allMetaData$lat, SR, pch=16, xaxt="s", las=1, ylab="", xlab="", col="darkgrey")
par(mar=c(4,0,1,1))
lines(x, p$fit, lwd=2)
lines(x, p$fit+p$se.fit*1.96, lwd=2, lty=2)
lines(x, p$fit-p$se.fit*1.96, lwd=2, lty=2)
mtext("Latitude (°)", 1, line=2.5)
mtext("Species Richness", 2, line=-1.5)
text(x=-25, y = 40, adj=4.5, labels="(a)", xpd=T, cex=1.3)

yhist = hist(SR, plot=FALSE, breaks=seq(from=min(SR), to=max(SR), length.out=15))
barplot(yhist$density, horiz=T, xaxt="n", col="white")

SpA = colSums(allSpeciesData[,-1]>0, na.rm=T)
par(mar=c(4,4,1,2), xaxt="s")
plot(log(sort(SpA, dec=T)), yaxt="n", xlab="", ylab="", las=1, pch=16)
axis(2, 0:5, signif(exp(0:5),1), las=1)
mtext("Abundance rank", 1, line=2.5)
mtext("Number of samples", 2, line=2.5)
text(x=0, y=5.4, adj=3.5, labels="(b)", xpd=T, cex=1.3)

# Simple map
plot(allMetaData$long, allMetaData$lat)

# Complex map
library(maptools)
library(raster)

brazil = readShapePoly("Mapdata/Brazil_biomes.shp")
names = c("Caatinga", "Cerrado", "Pantanal","Pampa","Amazônia", "Mata Atlântica")
map1 = readShapePoly("Mapdata/Political Map.shp")
LA = map1[map1$REGION=="Latin America",]
brazilRaster = raster(brazil)
res(brazilRaster) = 0.1
LAraster = raster(LA)
res(LAraster) = .1
brazil2 = rasterize(brazil, LAraster)

pdf("EUGCOMM_MAP.pdf", height=5, width=6.3, pointsize=10, family="Times")
par(oma=c(0,0,0,6))
plot(brazil2, legend=F, xlim=c(-95,-15), ylim=c(-36,20), col=topo.colors(10)[c(1,4,5,7,9,10)], box=F)
lines(LA)
points(allMetaData$long, allMetaData$lat, pch=16)

SpatialPolygonsRescale(layout.north.arrow(1), offset=c(-39,10), scale = 9, plot.grid=F)
arrows(-89.5, -34, -80.5, -34, code=3, length=0)
text(-85, -32.5, "1000 km")

plot(brazil2, legend.only=T, xlim=c(-93,-20), ylim=c(-36,25), col=topo.colors(10)[c(1,4,5,7,9,10)], 
     legend.args = list(axis=F, colors=1:6),
     axis.args = list(at=1:6, labels=names))
dev.off()


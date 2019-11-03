#######################################################################
############## - EUGCOMM Atlantic Forest HMSC analysis - ##############
#######################################################################

#### Model fitting ####

rm(list=ls())
gc()
memory.limit(64000)

#### Load packages ####
library(devtools)
library(withr)
#install_github("hmsc-r/HMSC", build_opts = c("--no-resave-data", "--no-manual"))
library(abind)
library(reshape2)
library(maptools)
library(corrplot)
library(ggplot2)
library(raster)
library(Hmsc)
library(fBasics)
library(fields)

setwd("F:/HelsinkiData23102019/euglossini")

# Load datafiles
Y = read.csv(file="data/Y.csv")
XData = read.csv(file="data/XData.csv")
TrData = read.csv(file="data/TrData.csv")
dfPi = read.csv(file="data/dfPi.csv")
xy = read.csv(file="data/xy.csv")
rownames(xy)=xy$X
xy=xy[,-1]

# Species synonyms
eco = which(colnames(Y)=="Euglossa_cordata")
Y$Euglossa_carolina = Y$Euglossa_carolina + Y$Euglossa_cordata
Y = Y[,-eco]
TrData = TrData[-eco,]

eci = which(colnames(Y)=="Eulaema_cingulata")
Y$Eulaema_marcii = Y$Eulaema_marcii + Y$Eulaema_cingulata
Y = Y[,-eci]
TrData = TrData[-eci]

eto = which(colnames(Y)=="Euglossa_townsendi")
Y$Euglossa_aratingae = Y$Euglossa_aratingae + Y$Euglossa_townsendi
Y = Y[,-eto]
TrData = TrData[-eto]

TrData = as.data.frame(TrData)
names(TrData) = "genus"

dim(Y)
dim(TrData)
dim(XData)
dim(dfPi)
dim(xy)

str(Y)
str(TrData)
str(XData)
str(dfPi)
str(as.matrix(xy))
head(xy)

# Set random levels
rL1 = HmscRandomLevel(sData = as.matrix(xy))
rL2 = HmscRandomLevel(units = unique(dfPi$SU))
rL1
rL2

# Set model formulae
XFormula = ~ method + effort + poly(altitude, 2, raw=TRUE) + MAT + MAP + Tseason + Pseason + 
  poly(forest., 2, raw=TRUE) + poly(lu_het, 2, raw=TRUE)

TrFormula= ~genus-1

#Setup HMSC models
m1 = Hmsc(Y = as.matrix((Y>0)*1), 
          XData = XData,  XFormula = XFormula,
          TrData = TrData, TrFormula = TrFormula,
          distr = "probit", 
          studyDesign = dfPi, 
          ranLevels = list(SA=rL1, SU=rL2))

Y2=Y
for(i in 1:nrow(Y2)){
  for(j in 1:ncol(Y2)){
    if(Y2[i,j]==0){
      Y2[i,j]=NA
    }}}

m2 = Hmsc(Y=as.matrix(log(Y2)), 
          XData = XData,  XFormula = XFormula,
          TrData = TrData, TrFormula = TrFormula,
          distr = "normal", 
          studyDesign = dfPi, 
          ranLevels = list(SA=rL1, SU=rL2))

head(m1$XData)
head(m2$XData)
cor(m1$XData[,2:9])

#### Sample MCMC ####
thin = 100
samples = 1000
adaptNf = ceiling(0.4*samples*thin)
transient = ceiling(0.5*samples*thin)
nChains = 2

a=Sys.time()
m1 = sampleMcmc(m1, samples = samples, thin = thin,
                adaptNf = rep(adaptNf, 2),
                transient = transient,
                nChains = nChains, nParallel = 1, updater=list(GammaEta=FALSE))
Sys.time()-a

save(m1,file="mMataPA_150K_2chains_Nov.Rdata")

a=Sys.time()
m2 = sampleMcmc(m2, samples = samples, thin = thin,
                adaptNf = rep(adaptNf, 2),
                transient = transient,
                nChains = nChains, nParallel = 2, updater=list(GammaEta=FALSE))
Sys.time()-a

#save(m2,file="mMataCond_150K_2chains_Sep.Rdata")

#### POSTPROCESSING ####

#### Presence-absence model ####
load("mMataPA_150K_2chains_Nov.Rdata")
m = m1

# Assess sampling performance
mpost = convertToCodaObject(m)

esBeta = effectiveSize(mpost$Beta)
summary(esBeta)

psrf = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
summary(psrf)

# Plot posterior trace
pdf("PAmod/alphaPost.pdf")
plot(mpost$Alpha[[1]])
dev.off()

pdf("PAmod/betaPost.pdf")
plot(mpost$Beta[,1:200])
dev.off()

# Explanatory power
predY = computePredictedValues(m, expected=T)
MF = evaluateModelFit(hM=m, predY = predY)

mean(MF$TjurR2, na.rm=T)
mean(MF$AUC, na.rm=T)

range(MF$TjurR2)
range(MF$AUC)

predYm = apply(predY, 1:2, mean)

plot(colSums(predYm,na.rm=T), colSums(m$Y))
lines(-1000:1000, -1000:1000)
cor(colSums(predYm), colSums(m$Y))^2

plot(rowMeans(predYm), rowMeans(m$Y))
lines(-1000:1000, -1000:1000)
cor(rowSums(predYm), rowSums(m$Y))^2

#### Variance partitioning ####
plotVP2=function (hM, VP, ...) 
{
  ng = dim(VP$vals)[1]
  leg = VP$groupnames
  for (r in 1:hM$nr) {
    leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
  }
  means = round(100 * rowMeans(VP$vals), 1)
  for (i in 1:ng) {
    leg[i] = paste(leg[i], " (mean = ", toString(means[i]), 
                   ")", sep = "")
  }
  mainTitle = substitute("")
  barplot(VP$vals, main = mainTitle, xlab = "", ylab = "Variance proportion", 
          las = 1, legend = leg, 
          ...)
}

cbind(m$covNames)
group = c(rep(1,3),2, rep(3,2) ,rep(4,4), rep(5,2), rep(6,2))
cbind(m$covNames, group)

groupnames = c("Baiting method", "Effort", "Altitude", "Climate", "Forest cover", "Landuse heterogeneity")
VP = computeVariancePartitioning(hM = m, group = group, groupnames = groupnames)
VP$R2T #Genus effects

VP$vals = VP$vals[,rev(order(colSums(VP$vals[1:6,])))]

pdf("PAmod/VarPartMata.pdf", height=3, width=7, family="Times")
par(mar=c(2,4,2,12))
plotVP2(m, VP, col=topo.colors(8), axisnames=F, args.legend=list(y=1.0, x=112, horiz=F, cex=0.75, bty="n"))
dev.off()

# Beta parameters
source("plotBeta.R")
library(fields)

pbeta=getPostEstimate(m,"Beta")

x11()
plotBeta(m,post=pbeta, param = "Support", covOrder="Vector", covVector=c(2:14), supportLevel = .95,
         SpeciesOrder = "Vector", SpVector = 1:58, spNamesNumbers=c(T,F), covNamesNumbers = c(T,F))

#### Gradient plots ####

x = 2
m$covNames
Gradient=constructGradient(m, focalVariable = "forest.", ngrid=50,
                           non.focalVariables = list(
                             method = list(1),
                             effort = list(1),
                             altitude = list(x),
                             MAT = list(x),
                             MAP = list(x),
                             Tseason = list(x),
                             Pseason = list(x),
                             lu_het = list(x)))
head(Gradient$XDataNew, 50)

predY = predict(m, Gradient=Gradient, expected=TRUE, predictEtaMean=FALSE)

plotGradient(m, Gradient, predY, measure="S", showData = T)
plotGradient(m, Gradient, predY, measure="Y", index=30, showData = T)

pdf("PAmod/ForestGradients.pdf", height=3.4, width=10, family="Times", pointsize=11)
par(mfrow=c(1,3), mar=c(4,4,2,1), oma=c(0,0,0,0))

xx = Gradient$XDataNew[,1]
cols = c(rgb(1, 0, 0, alpha = 0.5), rgb(0, 1, 0, alpha = 0.5), rgb(0, 0, 1, alpha = 0.5), rgb(0, 1, 1, alpha = 0.5))
cols = topo.colors(4, alpha=.5)

predS = lapply(predY, function(x) rowSums(x))
pred = apply(abind(predS, along = 2), c(1), quantile, 
             prob = c(.025, .5, .975), na.rm = TRUE)
plot(m$XData$forest., rowSums(m$Y), pch=16, lwd = 2, ylim=c(0, 30), col="grey", las=1,
     ylab="", xlab="")
mtext("Proportion forest cover", 1, line=2.5)
mtext("Species richness", 2, line=2.5, xpd=T)

lines(xx, pred[2, ], lwd=2)
polygon(c(xx, rev(xx)), c(pred[1, ], rev(pred[3, ])), 
        col = cols[3], border = FALSE)

text(x=0, y=30, adj=c(3.7), labels="(a)", xpd=T, cex=1.3)

# Per genus
predS = lapply(predY, function(x) rowSums(x[,which(m$TrData$genus=="Euglo")]))
egpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
plot(xx, egpred[2, ], type="l",lwd = 2, ylim=c(0, 20), las=1,
     ylab="", xlab="")
mtext("Proportion forest cover", 1, line=2.5)
polygon(c(xx, rev(xx)), c(egpred[1, ], rev(egpred[3, ])), 
        col = cols[1], border = FALSE)

predS = lapply(predY, function(x) rowSums(x[,which(m$TrData$genus=="Eulae")]))
elpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
lines(xx, elpred[2, ], lwd = 2)
polygon(c(xx, rev(xx)), c(elpred[1, ], rev(elpred[3, ])), 
        col = cols[2], border = FALSE)

predS = lapply(predY, function(x) rowSums(x[,which(m$TrData$genus=="Eufri")]))
efpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
lines(xx, efpred[2, ], lwd = 2)
polygon(c(xx, rev(xx)), c(efpred[1, ], rev(efpred[3, ])), 
        col = cols[3], border = FALSE)

predS = lapply(predY, function(x) rowSums(x[,which(m$TrData$genus=="Exaer")]))
expred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
lines(xx, expred[2, ], lwd = 2)
polygon(c(xx, rev(xx)), c(expred[1, ], rev(expred[3, ])), 
        col = cols[4], border = FALSE)

legendcols = c(rgb(1, 0, 0, alpha = 1),rgb(0, 1, 0, alpha = 1),rgb(0, 0, 1, alpha = 1),rgb(0, 1, 1, alpha = 1))
legendcols = topo.colors(4, alpha=0.75)
legend("topleft", c(expression(italic(Euglossa)), expression(italic(Eulaema)), expression(italic(Eufriesea)), expression(italic(Exaerete))),bty="n", cex=1.2 ,col=legendcols,lwd=2,lty=1)
mtext("Species richness", 2, line=2.5, xpd=T)

text(x=0, y=20, adj=c(3.7), labels="(b)", xpd=T, cex=1.3)

#Per species
cols = rev(divPalette(90, "Spectral"))[c(1:29, 62:90)]

predS = lapply(predY, function(x) x[,1])
egpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
plot(xx, egpred[2, ], type="l", lwd = 2, ylim=c(0,1), xlim=c(0,1), las=1,
    ylab="", xlab="", col=cols[1])

for(i in 2:ncol(m$Y)){
  predS = lapply(predY, function(x) x[,i])
  egpred = apply(abind(predS, along = 2), c(1), quantile, 
                 prob = c(.025, .5, .975), na.rm = TRUE)
  lines(xx, egpred[2, ], lwd = 2, col=cols[i])
}  

mtext("Proportion forest cover", 1, line=2.5)
mtext("Probability of occurrence", 2, line=2.5, xpd=T)

text(x=0, y=1, adj=c(3.7), labels="(c)", xpd=T, cex=1.3)

dev.off()


#### Abundance model ####
rm(list=ls())
load("mMataCond_150K_2chains_Sep.Rdata")
m = m2

#Assess sampling performance
mpost = convertToCodaObject(m)

esBeta = effectiveSize(mpost$Beta)
summary(esBeta)

psrf = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
summary(psrf)

esAlpha = effectiveSize(mpost$Alpha[[1]])
summary(esAlpha)

#Plot posterior trace
pdf("Amod/alphaPost.pdf")
plot(mpost$Alpha[[1]])
dev.off()

pdf("Amod/betaPost.pdf")
plot(mpost$Beta[,1:200])
dev.off()

#Explanatory power
predY = computePredictedValues(m, expected=T)
MF = evaluateModelFit(hM=m, predY = predY)

mean(MF$R2, na.rm=T)
range(MF$R2)

#### - Variance partitioning - ####
plotVP2=function (hM, VP, ...) 
{
  ng = dim(VP$vals)[1]
  leg = VP$groupnames
  for (r in 1:hM$nr) {
    leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
  }
  means = round(100 * rowMeans(VP$vals), 1)
  for (i in 1:ng) {
    leg[i] = paste(leg[i], " (mean = ", toString(means[i]), 
                   ")", sep = "")
  }
  mainTitle = substitute("")
  barplot(VP$vals, main = mainTitle, xlab = "", ylab = "Variance proportion", 
          las = 1, legend = leg, 
          ...)
}

cbind(m$covNames)
group = c(rep(1,3),2, rep(3,2) ,rep(4,4), rep(5,2), rep(6,2))
cbind(m$covNames, group)

groupnames = c("Baiting method", "Effort", "Altitude", "Climate", "Forest cover", "Landuse heterogeneity")
VP = computeVariancePartitioning(hM = m, group = group, groupnames = groupnames)

VP$R2T #Genus effects

VP$vals = VP$vals[,rev(order(colSums(VP$vals[1:6,])))]

tiff("Amod/VarPartMata.tiff", height=3, width=7, family="Times", units="in", res=400)
par(mar=c(2,4,2,12))
plotVP2(m, VP, col=topo.colors(8), axisnames=F, args.legend = list(y=1.0, x=118, horiz=F, cex=0.75, bty="n"))
dev.off()

#Beta parameters
source("plotBeta.R")

pbeta = getPostEstimate(m,"Beta")

x11()
plotBeta(m, post = pbeta, param = "Support", covOrder = "Vector", covVector=c(2:14), supportLevel = .75,
         SpeciesOrder = "Vector", SpVector = 1:58, spNamesNumbers=c(T,F), covNamesNumbers = c(T,F))

#### Gradient plots ####

x=2
m$covNames
Gradient = constructGradient(m, focalVariable = "forest.", ngrid=50,
                              non.focalVariables = list(
                              method=list(1),
                              effort=list(1),
                              altitude=list(x),
                              MAT=list(x),
                              MAP=list(x),
                              Tseason=list(x),
                              Pseason=list(x),
                              lu_het=list(x)))
head(Gradient$XDataNew, 50)

predY = predict(m, Gradient = Gradient, expected=TRUE, predictEtaMean=FALSE)

plotGradient(m, Gradient, predY, measure="S", showData = F)
plotGradient(m, Gradient, predY, measure="Y", index=30, showData = T)

tiff("Amod/ForestGradients.tiff", height=4, width=8, family="Times", units="in", res=400, pointsize=11)
par(mfrow=c(1,2), mar=c(4,4,2,1), oma=c(0,0,0,0))

xx = Gradient$XDataNew[,1]
cols = c(rgb(1, 0, 0, alpha = 0.5), rgb(0, 1, 0, alpha = 0.5), rgb(0, 0, 1, alpha = 0.5), rgb(0, 1, 1, alpha = 0.5))
cols = topo.colors(4, alpha=.5)

predS = lapply(predY, function(x) rowSums(x))
pred = apply(abind(predS, along = 2), c(1), quantile, 
             prob = c(.025, .5, .975), na.rm = TRUE)
plot(m$XData$forest.,rowSums(m$Y), pch=16, lwd = 2, ylim=c(40,120), col="grey", las=1,
     ylab="", xlab="")
mtext("Proportion forest cover", 1, line=2.5)
mtext("Abundance (sum log individuals)", 2, line=2.5, xpd=T)

lines(xx, pred[2, ], lwd=2)
polygon(c(xx, rev(xx)), c(pred[1, ], rev(pred[3, ])), 
        col = cols[3], border = FALSE)

text(x=0, y=120, adj=c(3.9), labels="(a)", xpd=T, cex=1.1)

# Per genus
predS = lapply(predY, function(x) rowSums(x[,which(m$TrData$genus=="Euglo")]))
egpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
plot(xx, egpred[2, ], type="l", lwd = 2, ylim=c(0, 120), las=1,
     ylab="", xlab="")
mtext("Proportion forest cover", 1, line=2.5)
polygon(c(xx, rev(xx)), c(egpred[1, ], rev(egpred[3, ])), 
        col = cols[1], border = FALSE)

predS = lapply(predY, function(x) rowSums(x[,which(m$TrData$genus=="Eulae")]))
elpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
lines(xx, elpred[2, ], lwd = 2)
polygon(c(xx, rev(xx)), c(elpred[1, ], rev(elpred[3, ])), 
        col = cols[2], border = FALSE)

predS = lapply(predY, function(x) rowSums(x[,which(m$TrData$genus=="Eufri")]))
efpred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
lines(xx, efpred[2, ], lwd = 2)
polygon(c(xx, rev(xx)), c(efpred[1, ], rev(efpred[3, ])), 
        col = cols[3], border = FALSE)

predS = lapply(predY, function(x) rowSums(x[,which(m$TrData$genus=="Exaer")]))
expred = apply(abind(predS, along = 2), c(1), quantile, 
               prob = c(.025, .5, .975), na.rm = TRUE)
lines(xx, expred[2, ], lwd = 2)
polygon(c(xx, rev(xx)), c(expred[1, ], rev(expred[3, ])), 
        col = cols[4], border = FALSE)

legendcols = c(rgb(1, 0, 0, alpha = 1), rgb(0, 1, 0, alpha = 1), rgb(0, 0, 1, alpha = 1), rgb(0, 1, 1, alpha = 1))
legendcols = topo.colors(4, alpha=0.75)
legend("topleft", c(expression(italic(Euglossa)), expression(italic(Eulaema)), expression(italic(Eufriesea)), expression(italic(Exaerete))), bty="n", cex=.8 , col=legendcols, lwd=2, lty=1)

text(x=0, y=120, adj=c(3.9), labels="(b)", xpd=T, cex=1.1)

dev.off()


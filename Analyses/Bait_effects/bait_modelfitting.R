##############################################################################
############## - EUGCOMM analysis of bait use and preferences - ##############
##############################################################################

rm(list=ls())
gc()
memory.limit(64000)

setwd("F:/HelsinkiData23102019/euglossini/baitmod")

library(devtools)
library(withr)
#install_github("hmsc-r/HMSC", build_opts = c("--no-resave-data", "--no-manual"))
library(Hmsc)
library(corrplot)
library(reshape2)

#### Model fitting

# Load datafiles
Y = read.csv(file="data/Y.csv")
XData = read.csv(file="data/XData.csv")
TrData = read.csv(file="data/TrData.csv")
rownames(TrData) = TrData$X
TrData = TrData[,-1]
dfPi = read.csv(file="data/dfPi.csv")

XData$Ska = as.factor(XData$Ska)
XData$Cin = as.factor(XData$Cin)
XData$Eug = as.factor(XData$Eug)
XData$MS = as.factor(XData$MS)
XData$MC = as.factor(XData$MC)
XData$BA = as.factor(XData$BA)
XData$VA = as.factor(XData$VA)

str(XData)
head(XData)
str(TrData)

#Set random effect structure
rL1 = HmscRandomLevel(units = unique(dfPi$SA))
rL2 = HmscRandomLevel(units = unique(dfPi$SU))

#Set model formulae
XFormula = ~ method + effort + Cin + Eug + MS + MC + BA + VA + Ska
TrFormula = ~ Cin + Eug + MS + MC + BA + VA + Ska

#Sample MCMC
thin = 200
samples = 1000
transient = .5*(thin*samples)
adaptNf = .4*(thin*samples)
nChains = 2

m = Hmsc(Y = as.matrix(Y>0)*1, 
         XData = XData,  XFormula = XFormula,
         TrData = TrData, TrFormula = TrFormula,
         distr = "probit", 
         studyDesign = dfPi, 
         ranLevels = list(SA=rL1, SU=rL2))
head(m$X)
head(m$Tr)
m$Y[1:5,1:5]
colnames(m$Y)==rownames(m$Tr)

a = Sys.time()
m = sampleMcmc(m, samples = samples, thin = thin,
             adaptNf = rep(adaptNf, m$nr),
             transient = transient,
             nChains = nChains, nParallel = 1, updater=list(GammaEta=FALSE))
Sys.time()-a

save(m, file="mBaits_Probit_200K_2chains_Nov.Rdata")

#### Load the model object ####
load(file="mBaits_Probit_200K_2chains_Nov.Rdata")

# Assess sampling performance
mpost = convertToCodaObject(m)

esBeta = effectiveSize(mpost$Beta)
summary(esBeta)

psrf = gelman.diag(mpost$Beta, multivariate=F)$psrf
summary(psrf)

esGamma = effectiveSize(mpost$Gamma)
summary(esGamma)

psrfGamma = gelman.diag(mpost$Gamma, multivariate=F)$psrf
summary(psrfGamma)

# Produce posterior trace plots
pdf("posterior_plots/betaPost.pdf")
plot(mpost$Beta[,1:200])
dev.off()

pdf("posterior_plots/gammaPost.pdf")
plot(mpost$Gamma)
dev.off()

# Explanatory power
predY = computePredictedValues(m, expected=T)
MF = evaluateModelFit(hM = m, predY = predY)
mean(MF$AUC, na.rm=T)
mean(MF$TjurR2, na.rm=T)
range(MF$AUC)
range(MF$TjurR2)

predYm = apply(predY, 1:2, mean)

plot(colSums(predYm,na.rm=T), colSums(m$Y))
lines(-1000:1000, -1000:1000)
cor(colSums(predYm), colSums(m$Y))^2

plot(rowMeans(predYm), rowMeans(m$Y))
lines(-1000:1000, -1000:1000)
cor(rowSums(predYm), rowSums(m$Y))^2

#### - Variance partitioning - ####
cbind(m$covNames)
group = c(rep(1,4), 2, rep(3,7))
cbind(m$covNames, group)
groupnames = c("Baiting method","Effort","Baits")

VP = computeVariancePartitioning(hM=m, group=group, groupnames=groupnames)

VP$R2T$Y #Variance in predicted occurrences explained by traits
VP$R2T$Beta #Variance in response to covariates explained by traits
mean(VP$R2T$Beta[6:12])

VP$vals=VP$vals[,rev(order(colSums(VP$vals[1:3,])))]
plotVariancePartitioning(hM=m, VP = VP)

# Beta parameters
pbeta = getPostEstimate(m, "Beta")
m$covNames
cbind(m$covNames[6:12], round(rowMeans(pbeta$mean[6:12,]),2))

summary(t(pbeta$mean))

tiff("baitBeta.tif", height=12, width=7, units="in", res=400, family="Times")
plotBeta(m, post=pbeta, param = "Mean", covOrder="Vector", covVector=c(6:12),
         SpeciesOrder = "Vector", SpVector=1:100, spNamesNumbers=c(T,F),
         covNamesNumbers = c(T,F), supportLevel = .95)
dev.off()

# Gamma parameters
pgamma = getPostEstimate(m, "Gamma")

plotGamma(m, post=pgamma, supportLevel=.75, param = "Mean",
          trOrder="Vector", trVector=c(2:8),
          covOrder="Vector", covVector=c(6:12), covNamesNumbers = c(T,T))

m$covNames
m$trNames
mat = pgamma$mean[6:12, 2:8]
smat = 2*pgamma$support[6:12, 2:8] - 1
supp = pgamma$support[6:12, 2:8]
supportLevel = .75
mat = mat * ((supp > supportLevel) + (supp < (1 - supportLevel)) > 0)

colnames(m$Tr)
names=c("Cineole", "Eugenol", "Methyl salicylate", "Methyl cinnamate", "Benzyl acetate", "Vanillin", "Skatole")
colnames(mat) = rownames(mat) = names

pdf("BaitPrefFig.pdf", height=5.5, width=11, family="Times")
par(mfrow=c(1,2), mar=c(5,7,3,0), xpd=F)
md = melt(pbeta$mean[6:12,])
plot(as.factor(md$Var1), md$value, xlab="", xaxt="n", ylab="Species response (probit scale)", las=1)
abline(h=0, lty=2, lwd=2)

axis(1, at=1:7, labels=rep("", 7))
text(x = seq(1,7, length.out = 7), par("usr")[3] - 
       +.5, srt = 45, adj = 1, cex = .8, labels = names, xpd = TRUE)
mtext("(a) Bait effects on species occurrence", 3, line=.5)

par(xpd=T)
corrplot(t(mat), is.cor=F,method="color", col=colorRampPalette(c("blue3", "white", "red3"))(200),
         mar=c(5,7,0,2), tl.cex=.6, tl.col="black", tl.pos="n",
         cl.align.text="r", cl.offset=0, cl.cex=.8, cl.length=11, cl.lim=c(-5,5), addgrid.col = "grey")
text(x = seq(1,7, length.out = 7), par("usr")[3] - 
       -1.5, srt = 45, adj = 1, cex = .8, labels = names, xpd = TRUE)
text(y = seq(1,7, length.out = 7), par("usr")[3] - 
       -1.5, srt = 0, adj = 1, cex = .8, labels = rev(names), xpd = TRUE)
mtext("Bait preference for", 2, line=5)
mtext("Response to inclusion of", 1, line=1)
mtext("(b) Effects of bait preference on bait response", 3, line=.5)

dev.off()

#### Gradient plots ####
m$covNames
type = 3
value = 1

Gradient = constructGradient(m, focalVariable = "MS",
                              non.focalVariables = list(
                              method=list(3, "Net"),
                              effort=list(1),
                              Ska=list(type, value),
                              Cin=list(type, value),
                              Eug=list(type, value),
                              MS=list(type, value),
                              MC=list(type, value),
                              BA=list(type, value),
                              VA=list(type, value)))

head(Gradient$XDataNew)

predY = predict(m, Gradient=Gradient, expected=TRUE, predictEtaMean=FALSE)

#Species richness
plotGradient(m, Gradient, predY, measure="S")

#Community-weighed mean trait value
m$trNames
plotGradient(m, Gradient, predY, measure = "T", index = 4)

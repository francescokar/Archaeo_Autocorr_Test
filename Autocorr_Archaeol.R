## Import dataset

library(GISTools)
setwd("F:/R-Code_Projects/Autocorrelation_Testing")
FB031<-readShapePoly("FB031_Surface")
RR_FB031<-readShapePoints("RR_FB031")

## Assess size and preservation proprtion of debris

max(RR_FB031$Dimensio_1)

tiff("Hist_Size_Debris.tif",width=480,height=360)
par(mar=c(4.5,4.5,3,1))
plot(hist(RR_FB031$Dimensio_1,breaks=30),col="#999999",xlim=c(0,800),xlab="size",
     ylab="frequency",main="Size of debris in FB031")
dev.off()

min(log(RR_FB031$Dimensio_1))
max(log(RR_FB031$Dimensio_1))

tiff("Hist_LogSize_Debris.tif",width=480,height=360)
plot(hist(log(RR_FB031$Dimensio_1),breaks=30),col="#666666",xlim=c(-2,7),
     xlab="log-size",ylab="frequency",main="Log-size of debris in FB031")
dev.off()


## Create plot

summary(RR_FB031$Dimensio_1)

tiff("RR_FB031.tif",width=600,height=480)
par(mar=c(0,0,4,0))
plot(FB031,lwd=1.8,col="#CCCCCC",main="FB031: intra-site debris distribution
     (log-size plotted)",cex.main=1.2)
plot(RR_FB031,cex=log(RR_FB031$Dimensio_1),pch=c(1,2)[RR_FB031$Preserva_1],
     lwd=1.5,col="red",add=T)
legend(-350,20,legend=c("Complete","Fragment"),pch=c(1,2),cex=1.2,col="red")
map.scale(-240,-110,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()


#### AUTOCORRELATION FOR CONTINUOUS VARIABLES ###

### Some useful references
### https://www.researchgate.net/post/Which_package_in_R_could_be_used_to_perform_Morans_I_test_for_spatial_autocorrelation


## GLOBAL METHODS ##

### MORAN'S I
## https://cran.r-project.org/web/packages/ape/vignettes/MoranI.pdf

## Moran.I function in library(ape)
# Non-spatial autocorrelation
library(ape)

# Compute distance weight (wij=dij, wii=1/dij)
RRcoord<-cbind(coordinates(RR_FB031))
RRdist<-dist(RRcoord,diag=T))
RRdist_mat<-as.matrix(RRdist)
RRweight<-1/RRdist_mat
diag(RRweight)<-0

# Moran's I, alternative weights
RRMoranI_1<-Moran.I(log(RR_FB031$Dimensio_1),RRdist_mat)
RRMoranI_2<-Moran.I(log(RR_FB031$Dimensio_1),RRweight)
write.table(capture.output(list("Wij=dij"=RRMoranI_1,"Wij=1/dij"=RRMoranI_2),
                           file="Moran.I_Results.txt"))

# correlogram.formula returns Moran's I for the different levels of taxonomic herarchy.
# difficult to adapt to archaeological spatial data

# to my knowledge, no option to compute Anselin's LISA


## moran.test, moran.mc and correlog functions in library spdep
library(spdep)

# Estimate nearest neighbour (k=4)
RRknn<-knearneigh(coordinates(RR_FB031),k=4,RANN=F)
RRnb<-knn2nb(RRknn)

# Moran's I test, with different spatal weight styles
RRmoran_W<-moran.test(log(RR_FB031$Dimensio_1),listw=nb2listw(RRnb,style="W"))
RRmoran_B<-moran.test(log(RR_FB031$Dimensio_1),listw=nb2listw(RRnb,style="B"))
write.table(capture.output(list(RRmoran_W,RRmoran_B)),file="moran.test_Results_B.txt")

# Moran's I test with Monte Carlo permutations & different spatial weight styles
RRmoran_mc_W<-moran.mc(log(RR_FB031$Dimensio_1),
                       listw=nb2listw(RRnb,style="W"),nsim=999)
RRmoran_mc_B<-moran.mc(log(RR_FB031$Dimensio_1),
                       listw=nb2listw(RRnb,style="B"),nsim=999)
write.table(capture.output(list(RRmoran_mc_W,RRmoran_mc_B)),
            file="moran.mc_Results.txt")

# Correlogram: Global Moran's I at multiple distances
RRcorrel<-sp.correlogram(neighbours=RRnb,var=log(RR_FB031$Dimensio_1),
                         order=5,method="I",style="W")

tiff("Correlogram_spdep.tif",width=600,height=480)
par(mar=c(5,5,5,3))
par(cex=1.3)
par(cex.main=1.2)
plot(RRcorrel,main="Moran's I correlogram log-size of RR")
dev.off()
# no direct option for plotting statistical significance

## GEARY'S C

# geary.test and sp.correlogram function in spdep package
library(spdep)
RRGeary<-geary.test(log(RR_FB031$Dimensio_1),nb2listw(RRnb))
write.table(capture.output(RRGeary),file="geary.test_Results.txt")

RRcorrel_C<-sp.correlogram(neighbours=RRnb,var=log(RR_FB031$Dimensio_1),
                           order=5,method="C",style="W")
tiff("Geary-correlogram_spdep.tif",width=600,height=480)
par(mar=c(5,5,5,3))
par(cex=1.3)
par(cex.main=1.2)
plot(RRcorrel_C,main="Geary's C correlogram log-size of RR")
dev.off()


## MORAN'S I, GEARY'S C AND CORRELOGRAMS WITH OTHER PACKAGES

# correlog function in pgrmess package
library(pgirmess)

RRcorrel<-correlog(coords=coordinates(RR_FB031),z=log(RR_FB031$Dimensio_1),
                       method="Moran",nbclass=5)
RRcorrel_alt<-correlog(coords=coordinates(RR_FB031),z=log(RR_FB031$Dimensio_1),
                       method="Geary",nbclass=5)
tiff("Correlog_pgirmess.tif",width=600,height=600)
par(mar=c(5,5,5,3))
par(mfrow=c(2,1))
par(cex=1.2)
plot(RRcorrel,main="Moran's I correlogram log-size of RR")
plot(RRcorrel_alt,main="Geary's I correlogram log-size of RR")
dev.off()

# USeful blog on correlograms
# https://www.r-bloggers.com/spatial-correlograms-in-r-a-mini-overview/

## gearymoran function in ade4 packages
library(ade4)

# Different weights tested, as in Moran.I
RR_gearymoran_1<-gearymoran(RRdist_mat,log(RR_FB031$Dimensio_1),nrepet=999)
RR_gearymoran_2<-gearymoran(RRweight,log(RR_FB031$Dimensio_1),nrepet=999)
write.table(capture.output(list("Wij=dij"=RR_gearymoran_1,
                                "Wij=1/dij"=RR_gearymoran_2)),file="gearymoran_Results.txt")


## MANTEL TEST AND CORRELOGRAM

# sp.mantel.mc function in spdep package (Mantel-Hubert Spatial General Cross Product Statistic)
library(spdep)
RR_spmantel<-sp.mantel.mc(log(RR_FB031$Dimensio_1),listw=nb2listw(RRnb,style="W"),nsim=999,type="moran")
write.table(capture.output(RR_spmantel),file="sp.mantel.mc_Results_W.txt")
plot.mc.sim(RR_spmantel)

# mantel.test function in ape package
library(ape)

RRsize<-dist(log(RR_FB031$Dimensio_1),diag=T)
RRsize_mat<-as.matrix(RRsize)

# Mantel test
RRMantel<-mantel.test(RRdist_mat,RRsize_mat,nperm=999)
write.table(capture.output(RRMantel),file="mantel.test_Results.txt")


# mantel.correlog function in vegan package

library(vegan)
RRMantel_correl<-mantel.correlog(RRsize,XY=coordinates(RR_FB031),
                                 n.class=5,nperm=999,cutoff=300)
tiff("Mantel_correlog.tif",width=600,height=480)
par(mar=c(5,5,5,3))
par(cex=1.2)
plot(RRMantel_correl,alpha=0.05)
title("Mantel correlogram log-size of RR")
dev.off()


## LOCAL METHODS ##

# local Moran's I and Getis-Ord Gi in spdep (localmoran & localG)
library(spdep)
RRlm<-localmoran(log(RR_FB031$Dimensio_1),listw=nb2listw(RRnb,style="W"))
RRlm_I1<-(RRlm[,1]<=-1)+0
RRlm_I2<-(RRlm[,1]>-1 & RRlm[,1]<=0)*2
RRlm_I3<-(RRlm[,1]>0 & RRlm[,1]<=1)*3
RRlm_I4<-(RRlm[,1]>1)*4
RRlm_Iclass<-RRlm_I1+RRlm_I2+RRlm_I3+RRlm_I4
RRlm_pch<-c(15,17,18,19)
RRlm_pr<-(RRlm[,5]<=0.05)+0
RRlm_pr2<-(RRlm[,5]>0.05)*2
RRlm_Prclass<-RRlm_pr+RRlm_pr2
RRlm_col<-c("#31A354","#DE2D26")

tiff("Localmoran.tif",width=600,height=480)
par(mar=c(0,2,5,0))
plot(FB031,lwd=1.9,col="#CCCCCC",main="FB031: Local-Moran's I",cex.main=1.5)
plot(RR_FB031,cex=1.5,pch=RRlm_pch[RRlm_Iclass],col=RRlm_col[RRlm_Prclass],add=T)
legend(-350,20,legend=c("Pr < 0.05","Pr > 0.05"),pch=15,pt.cex=1.2,
       col=c("#31A354","#DE2D26"))
legend(-350,140,legend=c("I < -1","I > -1 < 0","I > 0 < 1","I > 1"),pt.cex=1.2,
       pch=RRlm_pch)
map.scale(-240,-110,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()

RRlG<-localG(log(RR_FB031$Dimensio_1),listw=nb2listw(RRnb,style="W"))
RRlG_1<-(RRlG<=-1)+0
RRlG_2<-(RRlG>-1 & RRlG<=0)*2
RRlG_3<-(RRlG>0 & RRlG<=1)*3
RRlG_4<-(RRlG>1)*4
RRlG_class<-RRlG_1+RRlG_2+RRlG_3+RRlG_4
RRlG_pch<-c(15,17,18,19)
tiff("LocalG.tif",width=600,height=480)
par(mar=c(0,2,5,0))
plot(FB031,lwd=1.8,col="#CCCCCC",main="FB031: Getis-Ord Gi",cex.main=1.5)
plot(RR_FB031,cex=1.5,pch=RRlG_pch[RRlG_class],add=T)
legend(-350,50,legend=c("Gi < -1","Gi > -1 < 0","Gi > 0 < 1","Gi > 1"),pt.cex=1.2,
       pch=RRlG_pch)
map.scale(-240,-110,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()


## SEMIVARIOGRAM ##

# variogram function in gstat package
library(gstat)
RRvariogram<-variogram(log(RR_FB031$Dimensio_1)~1,RR_FB031,cutoff=150,width=15)
tiff("Variogram_gstat.tif",width=480,height=480)
par(mar=c(5,5,3,2))
plot(gamma~dist,RRvariogram,pch=16,xlab="distance(cm)",ylab="semivariance",
     col="black",main="Variogram log-size of RR")
lines(RRvariogram$dist,RRvariogram$gamma,lty=2,lwd=1.5)
dev.off()

RRdirectvar<-variogram(RR_FB031$Dimensio_1~1,RR_FB031,alpha=c(0,45,90,135),
                       cutoff=150,width=15)
tiff("Directvariogram_gstat.tif",width=480,height=480)
plot(RRdirectvar,as.table=T,pch=16,lty=2,lwd=1.5,xlab="distance(cm)",
     main="Directional variogram log-size of RR",col="black")
dev.off()

# variog function in geoR package
library(geoR)
RRvariog<-variog(coords=coordinates(RR_FB031),data=log(RR_FB031$Dimensio_1),max.dist=150)
tiff("Variog_geoR.tif",width=480,height=480)
par(mar=c(5,5,3,2))
plot(RRvariog,main="Omnidirectional variogram of log-size of RR",pch=20,
     xlab="distance(cm)",ylab="semivariance")
lines(RRvariog,lty=2,lwd=1.5)
dev.off()

# directional variogram also available 
# (option direction=[0,180],unit.angle="degrees")


### AUTOCORRELATION FOR CATEGORICAL VARIABLES ####

## CROSS-L FUNCTION ANALYSIS ##

# envelope function in spatstat package
library(spatstat)
FB031owin<-as.owin(FB031)
RRppp<-as.ppp(coordinates(RR_FB031),FB031owin)
marks(RRppp)<-RR_FB031$Preserva_1
plot(RRppp)
RRLcross<-envelope(RRppp,fun=Lcross,i="Complete",j="Fragment",
                   nsim=999,correction="Ripley")
tiff("Lcross.tif",width=480,height=480)
par(mar=c(5,5,5,3))
plot(RRLcross,main="Cross-L Function: Complete, Fragment")
dev.off()

RR_dclftest<-dclf.test(RRppp,Lcross,i="Complete",j="Fragment",
                       rinterval=c(5,20),nsim=999)
write.table(capture.output(RR_dclftest),file="dclf.test_Result.txt")
# Works only with points!

## JOIN-COUNT (BW) STATISTICS

# joincount.test, joincount.mc & joincount.multi functions in spdep package
library(spdep)
RRjc<-joincount.test(RR_FB031$Preserva_1,listw=nb2listw(RRnb))
write.table(capture.output(RRjc),file="joincount.test_Results.txt")
# same colour and different colour statistics
RRjc_multi<-joincount.multi(RR_FB031$Preserva_1,listw=nb2listw(RRnb))
write.table(capture.output(RRjc_multi),file="joincount.multi_Results.txt")
# permutation test j.c. statistics
RRjc_mc<-joincount.mc(RR_FB031$Preserva_1,listw=nb2listw(RRnb),nsim=999)
write.table(capture.output(RRjc_mc),file="joincount.mc_Results.txt")


## LOCAL INDICATOR FOR CATEGORICAL DATA ##

# Provisional example

RR_FB031_categ<-as.numeric(RR_FB031$Preserva_1)-1

for(i in 1:length(RR_FB031_categ)){
  xi<-(RR_FB031_categ[1:i])
  xj<-RR_FB031_categ
  w<-(RRweight[,1:i])
  bb<-w*xj
  ww<-w*(1-xj)}

sbb<-colSums(bb)
BB_RRLICD<-RR_FB031_categ*s
sww<-colSums(ww)
WW_RRLICD<-(1-RR_FB031_categ)*sww
summary(BB_RRLICD)
summary(WW_RRLICD)

# plot the LICD results
RR_BB<-(BB_RRLICD>=0.3029)*3
RR_WW<-(WW_RRLICD>=0.11200)*2
RR_BW<-(BB_RRLICD<0.3029 & WW_RRLICD<0.11200)+0
RR_BBWW<-RR_BB+RR_WW+RR_BW
RRJCS_pch<-c(4,16,1)
tiff("LICD.tif",width=600,height=480)
par(mar=c(0,2,5,0))
plot(FB031,lwd=1.8,col="#CCCCCC",main="FB031: Local-JCS",cex.main=1.5)
plot(RR_FB031,cex=1.5,pch=RRJCS_pch[RR_BBWW],add=T)
legend(-350,10,legend=c("Not significant","BB","WW"),pt.cex=1.2,
       pch=RRJCS_pch)
map.scale(-240,-110,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()

# compare with Complete/Fragment distribution
tiff("Preservation.tif",width=600,height=480)
par(mar=c(0,2,5,0))
plot(FB031,lwd=1.8,col="#CCCCCC",main="FB031: Preservation",cex.main=1.5)
plot(RR_FB031,cex=1.5,pch=c(16,1)[RR_FB031$Preserva_1],add=T)
legend(-350,10,legend=c("Complete","Fragment"),pt.cex=1.2,
       pch=RRJCS_pch)
map.scale(-240,-110,len=200,units="meters",ndivs=2,subdiv=1)
dev.off()
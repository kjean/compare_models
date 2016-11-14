##### 
## print output with FOI effect


rm(list=ls(all=TRUE))


library(maptools)
#library(sp) 
library(shapefiles)
library(fields)
library(Hmisc)
library(gplots)
library(ggplot2)
library(RColorBrewer)

# 0: labo, 1: portable
computer<-3
if (computer==0) { homedir<-"Y:/Kevin/"
} else if (computer==1) { homedir<-"C:"
} else if (computer==2) { homedir<-"//fi--didenas1/yf/Kevin/"
} else if (computer==3) { homedir<-"C:/Users/Kevin JEAN/Desktop/Recherche/Yellow_Fever/Kevin_YF_DIDE_share_drive/"
} else if (computer==4) { homedir<-"E:/Imperial/YF share drive/Kevin/"
} else if (computer==5) { homedir<-"/media/kevinNFS/"
}

shpdir = paste(homedir,"shapefiles/gadm2/",sep="")
commondir = paste(homedir,"re-fit_models_2016/", sep="")
currdir = paste(homedir,"re-fit_models_2016/compare_models", sep="")



shp0 = readShapePoly(paste(shpdir, "Africa_adm0.shp",sep=""))
shp1 = readShapePoly(paste(shpdir, "Africa_adm1.shp",sep=""))

shp1$adm0_adm1<-paste(shp1$ISO, shp1$ID_1, sep="_")
shp1 = shp1[order(shp1$adm0_adm1),]

setwd(currdir)
#load("updated_model_without_herd_data.Rdata")


########
# FOI
foi_dir =paste(commondir,'YF_FOI_burden_model/MCMC_FOI_nb_runs=1000_20161025', sep='')
setwd(foi_dir)

load("../workspace_FOI_model_mcmc.Rdata")

vec = 2:250
foi_africa= NULL
file=paste('FOI_bb',1,'.csv', sep='' )
foi_africa = read.csv(file, stringsAsFactors=F)
foi_africa = foi_africa[,-1]
dim(foi_africa)


for(i in vec) {
  file=paste('FOI_bb',i,'.csv', sep='' )
  if (file.exists(file)){
    foi_africa_tmp = read.csv(file, stringsAsFactors=F)
    foi_africa_tmp = foi_africa_tmp[,-1]
    foi_africa = cbind(foi_africa, foi_africa_tmp)
  }
}
rm(foi_africa_tmp)
dim(foi_africa)

foi_med = apply(as.matrix(foi_africa),1, quantile, probs=0.5)
foi_low = apply(as.matrix(foi_africa),1, quantile, probs=0.025)
foi_sup = apply(as.matrix(foi_africa),1, quantile, probs=0.975)
foi_quantiles = data.frame(median =foi_med , low=foi_low, sup=foi_sup)





#####
#R0
R0_dir = paste(commondir,"YF_R0_burden_model/finalMCMC_R0_nb_runs=1000_20161028",sep="")
setwd(R0_dir)

R0_africa= NULL
file=paste('R0_bb',1,'.csv', sep='' )
R0_africa = read.csv(file, stringsAsFactors=F)
R0_africa = R0_africa[,-1]
dim(R0_africa)


for(i in vec) {
  file=paste('R0_bb',i,'.csv', sep='' )
  if (file.exists(file)){
    R0_africa_tmp = read.csv(file, stringsAsFactors=F)
    R0_africa_tmp = R0_africa_tmp[,-1]
    R0_africa = cbind(R0_africa, R0_africa_tmp)
  }
}
rm(R0_africa_tmp)
dim(R0_africa)

R0_med = apply(as.matrix(R0_africa),1, quantile, probs=0.5)
R0_low = apply(as.matrix(R0_africa),1, quantile, probs=0.025)
R0_sup = apply(as.matrix(R0_africa),1, quantile, probs=0.975)
R0_quantiles = data.frame(median =R0_med , low=R0_low, sup=R0_sup)





######
## plotting the results



setwd(paste(currdir, "/outputs_figures/", sep=""))

png("map_FOI_R0.png", width=18,height=9,units="in",res=400)

par(mfrow=c(1,2))


range(foi_med)
foi_percent = 100*foi_med
range(foi_percent)
log_foi_per = log10(foi_percent)
range(log_foi_per)
ticks = c(1e-4, 1e-3, 1e-2, 1e-1,1,10)

mybreaks= seq(min(log_foi_per)-0.000001, max(log_foi_per)+0.000001, length.out=101)
#mycols =  tim.colors(length(mybreaks)-1)
mycols = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)) 
#mycols = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) 
mm = match(shp1$adm0_adm1,dat$adm0_adm1)
vcols = findInterval(log_foi_per,mybreaks)

par(oma=c(4,4,4,4))

par(mar=c(2,2,4,4))
#png("map_FOI.png", width=16,height=12,units="in",res=200)
plot(shp0, xlim=c(-15,50),ylim=c(-20,30), main = "A. Estimated FOI (%)", cex.main=1.7)
mm0 = match(shp0$ISO,c34) #
plot(shp0[is.na(mm0),],col="grey80",add=T) 
plot(shp1[!is.na(mm),],col=mycols[vcols], xlim=c(-15,45),ylim=c(-20,30) , lty=0, add=T)#,border=mycols[vcols])
plot(shp0, lwd=2, add=T)
image.plot(legend.only=T,breaks=mybreaks,col=mycols,zlim=c(0,1), axis.args=list( at=log10(ticks), labels=ticks), legend.cex=3)
#dev.off()


## try in log R0
R0_med_20 = ifelse(R0_med>20, 20,R0_med)
R0_med_10 = ifelse(R0_med>10, 10,R0_med)
log_R0 = log(R0_med_20)
log_log_R0 = log(log(R0_med_10))


mybreaks2= seq(min(log_log_R0)-0.00001, log(log(10))+000001, length.out=101)
#mycols2 =  heat.colors(length(mybreaks2)-1)
mycols2 = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)) 
#mycols2 = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) 
mm = match(shp1$adm0_adm1,dat$adm0_adm1)
vcols2 = findInterval(log_log_R0,mybreaks2)
#vcols2=101-vcols2



ticks2 = c(1.0001, 1.001,1.01,  1.1,  1.5, 10)
labticks2 = c(1.0001,1.001,1.01,  1.1,  1.5, "10+")


#par(oma=c(4,2,4,4))
par(mar=c(2,2,4,4))
#par(oma=c(2,4,4,4))
#par(mar=c(2,4,4,6))
#png("map_log_R0.png", width=16,height=12,units="in",res=200)
plot(shp0, xlim=c(-15,50),ylim=c(-20,30),main = "B. Estimated R0", cex.main=1.7)
mm0 = match(shp0$ISO,c34) #
plot(shp0[is.na(mm0),],col="grey80",add=T) 
plot(shp1[!is.na(mm),],col=mycols2[vcols2], xlim=c(-15,45),ylim=c(-18,30) , lty=0, add=T)#,border=mycols[vcols])
plot(shp0, lwd=2, add=T)
#par( lab.cex=2.0)
image.plot(legend.only=T,breaks=mybreaks2,col=mycols2,zlim=c(0,1), 
           axis.args=list( at=log(log(ticks2)), labels=labticks2), legend.cex=3)
dev.off()





#############################################################################
########### R0 translated into FOI
#############################################################################

currdir = paste(homedir,"script_MCMC_herd_immunity/last_version_11_Feb",sep="")
setwd(paste(currdir,'/Thu_Feb_11_20164compartiments_finalMCMC_nb_runs=1000', sep=''))

load("../herd_data_foi_const_new_surv.Rdata")
source(paste(homedir,"script_MCMC_herd_immunity/functions_MCMC_R0_4_compartiments.R", sep=""))



R0_africa= NULL
file=paste('R0_bb',1,'.csv', sep='' )
R0_africa = read.csv(file, stringsAsFactors=F)
R0_africa = R0_africa[,-1]
dim(R0_africa)


for(i in vec) {
  file=paste('R0_bb',i,'.csv', sep='' )
  if (file.exists(file)){
    R0_africa_tmp = read.csv(file, stringsAsFactors=F)
    R0_africa_tmp = R0_africa_tmp[,-1]
    R0_africa = cbind(R0_africa, R0_africa_tmp)
  }
}
rm(R0_africa_tmp)
dim(R0_africa)

R0_med = apply(as.matrix(R0_africa),1, quantile, probs=0.5)
R0_low = apply(as.matrix(R0_africa),1, quantile, probs=0.025)
R0_sup = apply(as.matrix(R0_africa),1, quantile, probs=0.975)
R0_quantiles = data.frame(median =R0_med , low=R0_low, sup=R0_sup)

dim(pop.moments.whole)
pop.moments.whole = data.frame(adm = dn1, pop.moments.whole)
foi_med_transform = foi_whole_prevac(adm= dn1, R0 =R0_med,  pop.moments=pop.moments.whole, polydeg=6)
range(foi_med_transform)




##### plot foi from FOI and R0 models



setwd(paste(homedir, "MCMC_compare_both_models/", sep=""))



png("map_FOI_FOI_from_R0.png", width=16,height=9,units="in",res=100)

par(mfrow=c(1,2))


range(foi_med)
foi_percent = 100*foi_med
range(foi_percent)
log_foi_per = log10(foi_percent)
range(log_foi_per)
ticks = c(1e-4, 1e-3, 1e-2, 1e-1,1,10)

mybreaks= seq(min(log_foi_per)-0.000001, max(log_foi_per)+0.000001, length.out=101)
#mycols =  tim.colors(length(mybreaks)-1)
mycols = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)) 
#mycols = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) 
mm = match(shp1$adm0_adm1,dat$adm0_adm1)
vcols = findInterval(log_foi_per,mybreaks)

par(oma=c(4,4,4,2))
par(mar=c(2,4,4,4))
#png("map_FOI.png", width=16,height=12,units="in",res=200)
plot(shp0, xlim=c(-15,45),ylim=c(-20,30), cex.main=1.7)
#text(-14, 35, "A.", cex = 3, font = 2)
mm0 = match(shp0$ISO,c34) #
plot(shp0[is.na(mm0),],col="grey80",add=T) 
plot(shp1[!is.na(mm),],col=mycols[vcols], xlim=c(-15,45),ylim=c(-20,30) , lty=0, add=T)#,border=mycols[vcols])
plot(shp0, lwd=2, add=T)
image.plot(legend.only=T,breaks=mybreaks,col=mycols,zlim=c(0,1), axis.args=list( at=log10(ticks), labels=ticks), legend.cex=3)
#dev.off()






range(foi_med_transform)
foi_percent_transform = 100*foi_med_transform
range(foi_percent_transform)
log_foi_per_transform = log10(foi_percent_transform)
range(log_foi_per_transform)
ticks2 = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1,1,10)


log_foi_per_transform = ifelse (log_foi_per_transform<min(log_foi_per), min(log_foi_per), log_foi_per_transform)

mybreaks= seq(min(log_foi_per_transform)-0.000001, max(log_foi_per)+0.000001, length.out=101)
#mycols =  tim.colors(length(mybreaks)-1)
mycols = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)) 
#mycols = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) 
mm = match(shp1$adm0_adm1,dat$adm0_adm1)
vcols = findInterval(log_foi_per_transform,mybreaks)

#par(oma=c(4,4,4,2))
par(mar=c(2,4,4,6))
#png("map_FOI.png", width=16,height=12,units="in",res=200)
plot(shp0, xlim=c(-15,45),ylim=c(-20,30), cex.main=1.7)
#text(-14, 36, "B.", cex = 3, font = 2)
mm0 = match(shp0$ISO,c34) #
plot(shp0[is.na(mm0),],col="grey80",add=T) 
plot(shp1[!is.na(mm),],col=mycols[vcols], xlim=c(-15,45),ylim=c(-20,30) , lty=0, add=T)#,border=mycols[vcols])
plot(shp0, lwd=2, add=T)
image.plot(legend.only=T,breaks=mybreaks,col=mycols,zlim=c(0,1), axis.args=list( at=log10(ticks2), labels=ticks2), legend.cex=3)

dev.off()





#############################################################################
########### R0 translated into FOI
#############################################################################



###### plot R0 and CVC
load("Y:/Kevin/script_MCMC/vc3d.Rdata")
load("Y:/Kevin/script_MCMC/p_prop3d.Rdata")
load("Y:/Kevin/script_MCMC/P_tot_2d.Rdata")
dim(vc3d)

vac_cov_2014 = rep(NA, dim(vc3d)[1])
for (adm in 1:dim(vc3d)[1]){
  vac_cov_2014[adm] = sum(vc3d[adm, 2014-1939,]*p_prop3d[adm, 2014-1939,], na.rm=T)
}


R0_med_20 = ifelse(R0_med>20, 20,R0_med)
R0_med_10 = ifelse(R0_med>10, 10,R0_med)
log_R0 = log(R0_med_20)
log_log_R0 = log(log(R0_med_10))


cvc = 1-1/R0_med
range(cvc)
cvc_reached = rep(NA, dim(vc3d)[1])
cvc_reached[vac_cov_2014>cvc]=T
cvc_reached[vac_cov_2014<=cvc]=F
table(cvc_reached)
dimnames(vc3d)
names(cvc_reached) = dimnames(vc3d)[[1]]



mybreaks2= seq(min(log_log_R0)-0.00001, log(log(10))+000001, length.out=101)
#mycols2 =  heat.colors(length(mybreaks2)-1)
mycols2 = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)) 
#mycols2 = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) 
mm = match(shp1$adm0_adm1,dat$adm0_adm1)
vcols2 = findInterval(log_log_R0,mybreaks2)
#vcols2=101-vcols2



ticks2 = c(1.01, 1.05, 1.1, 1.25, 2, 10)
labticks2 = c(1.01, 1.05, 1.1, 1.25, 2, "10+")



setwd(paste(homedir, "MCMC_compare_both_models/", sep=""))

png("map_R0_and_CVC.png", width=16,height=9,units="in",res=100)
par(mfrow = c(1,2))
#par(oma=c(2,4,4,4))
par(mar=c(2,4,4,6))
#png("map_log_R0.png", width=16,height=12,units="in",res=200)
plot(shp0, xlim=c(-15,45),ylim=c(-20,30))
mm0 = match(shp0$ISO,c34) #
plot(shp0[is.na(mm0),],col="grey80",add=T) 
plot(shp1[!is.na(mm),],col=mycols2[vcols2], xlim=c(-15,45),ylim=c(-18,30) , lty=0, add=T)#,border=mycols[vcols])
plot(shp0, lwd=2, add=T)
#par( lab.cex=2.0)
image.plot(legend.only=T,breaks=mybreaks2,col=mycols2,zlim=c(0,1), 
           axis.args=list( at=log(log(ticks2)), labels=labticks2), legend.cex=3)


plot(shp0, xlim=c(-15,45),ylim=c(-20,30))
mm0 = match(shp0$ISO,c34) #
plot(shp0[is.na(mm0),],col="grey80",add=T) 
plot(shp1[!is.na(mm),],col=mycols2[vcols2], xlim=c(-15,45),ylim=c(-18,30) , lty=0, add=T)#,border=mycols[vcols])
mm_cvc = match(names(cvc_reached)[cvc_reached==T], shp1$adm0_adm1)
plot(shp1[mm_cvc,], col="#636363", add=T)
plot(shp0, lwd=2, add=T)
#par( lab.cex=2.0)
image.plot(legend.only=T,breaks=mybreaks2,col=mycols2,zlim=c(0,1), 
           axis.args=list( at=log(log(ticks2)), labels=labticks2), legend.cex=3)
legend("bottomleft", fill = "#636363", legend="CVC atteinte en 2014", cex = 1.5)
dev.off()

library(maptools)
library(sp) 
library(shapefiles)
library(fields)
library(Hmisc)
library(fields)
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
}

shpdir = paste(homedir,"shapefiles/gadm2/",sep="")
commondir = paste(homedir,"re-fit models 2016/", sep="")
currdir = paste(homedir,"re-fit_models_2016/", "compare_models/", sep="")
outdir = paste0(currdir, "outputs_figures/")

setwd(currdir)



####### ####### ####### ####### ####### 
####### load shapefiles


# shp0 = readShapePoly(paste(shpdir, "Africa_adm0.shp",sep=""))
# shp1 = readShapePoly(paste(shpdir, "Africa_adm1.shp",sep=""))
# 
# shp1$adm0_adm1<-paste(shp1$ISO, shp1$ID_1, sep="_")
# shp1 = shp1[order(shp1$adm0_adm1),]

c34 = c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ERI", "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GNQ", "KEN", "LBR", "MLI", "MRT", "NER", "NGA", "RWA", "SDN", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZMB")
country34 = c("Angola","Burundi","Benin","Burkina Faso","Central African Republic","Côte d'Ivoire","Cameroon","Democratic Republic of the Congo","Republic of Congo","Eritrea","Ethiopia","Gabon","Ghana","Guinea","Gambia","Guinea-Bissau","Equatorial Guinea","Kenya","Liberia","Mali","Mauritania","Niger","Nigeria","Rwanda","Sudan","Senegal","Sierra Leone","Somalia", "South Sudan", "Chad","Togo","Tanzania","Uganda", "Zambia") 



####### CFR distrib
set.seed(101)
prop.death.all.cases = rlnorm(1000, meanlog= -3.00, sdlog=0.442)
prop.death.all.cases = ifelse(prop.death.all.cases>1,1,prop.death.all.cases)
prop.death.all.cases = rep(prop.death.all.cases, each = 34)

# severe cases distrib
set.seed(101)
prop.severe = rlnorm(1000,meanlog=-2.222, sdlog=0.427)
prop.severe = ifelse(prop.severe>1,1, prop.severe)
prop.severe = rep(prop.severe, each = 34)


#######################################################
# import burden FOI model

foi_dir = paste(homedir, "script_MCMC/post_GAVI_runs_FOI_model/best_estimates/", sep="")
tab = read.csv(paste0(foi_dir,"cases_by_year_adm0_nb_runs=1000.csv"), h=T)
dim(tab)
colnames(tab)
head(tab[,1:4])
tab = tab[,-2]
colnames(tab)
tab[,2:ncol(tab)]=prop.severe*tab[,2:ncol(tab)] # apply CFR


burd_mean_foi = burd_med_foi = burd_inf_foi = burd_sup_foi = NULL
for(adm in c34){
  print(adm)
  temp = tab[tab$adm0 == adm,]
  mean_tmp = apply(temp[,-1], 2, mean)
  med_tmp = apply(temp[,-1], 2, median)
  low_tmp = apply(temp[,-1], 2, quantile, probs=0.025)
  sup_tmp = apply(temp[,-1], 2,  quantile, probs=0.975)
  burd_mean_foi = rbind(burd_mean_foi, mean_tmp)
  burd_med_foi = rbind(burd_med_foi, med_tmp)
  burd_inf_foi = rbind(burd_inf_foi, low_tmp)
  burd_sup_foi = rbind(burd_sup_foi, sup_tmp)
}



#######################################################
# import 1984-2013 burden R0 model
r0_dir = paste(homedir, "script_MCMC_herd_immunity/post_GAVI_round/Thu_Mar_10_20164compartiments_finalMCMC_nb_runs=1000/best_estimates/", sep="")
tab = read.csv(paste0(r0_dir,"cases_by_year_adm0_nb_runs=1000.csv"), h=T)
dim(tab)
colnames(tab)
tab= tab[,-2]
tab[,2:ncol(tab)]=prop.severe*tab[,2:ncol(tab)] # apply CFR

burd_mean_R0 = burd_med_R0 = burd_inf_R0 = burd_sup_R0 = NULL
for(adm in c34){
  print(adm)
  temp = tab[tab$adm0 == adm,]
  mean_tmp = apply(temp[,-1], 2, mean)
  med_tmp = apply(temp[,-1], 2, median)
  low_tmp = apply(temp[,-1], 2, quantile, probs=0.025)
  sup_tmp = apply(temp[,-1], 2,  quantile, probs=0.975)
  burd_mean_R0 = rbind(burd_mean_R0, mean_tmp)
  burd_med_R0 = rbind(burd_med_R0, med_tmp)
  burd_inf_R0 = rbind(burd_inf_R0, low_tmp)
  burd_sup_R0 = rbind(burd_sup_R0, sup_tmp)
}



######################################################
### select burden years

### select 1984-2013
colnames(burd_mean_foi)[1:30]
burd_foi = burd_mean_foi[,1:30]
burd_foi_inf = burd_inf_foi[,1:30]
burd_foi_sup = burd_sup_foi[,1:30]

burd_R0 = burd_mean_R0[,1:30]
burd_R0_inf = burd_inf_R0[,1:30]
burd_R0_sup = burd_sup_R0[,1:30]



burd_tot_foi = colSums(burd_foi)
burd_tot_foi_inf= colSums(burd_foi_inf)
burd_tot_foi_sup = colSums(burd_foi_sup)

burd_tot_R0 = colSums(burd_R0)
burd_tot_R0_inf= colSums(burd_R0_inf)
burd_tot_R0_sup = colSums(burd_R0_sup)

burd_tot = rbind(burd_tot_foi, burd_tot_R0)
colnames(burd_tot) = 1984:2013
burd_tot_inf = rbind(burd_tot_foi_inf, burd_tot_R0_inf)
burd_tot_sup = rbind(burd_tot_foi_sup, burd_tot_R0_sup)



png(paste0(outdir,"compare_burden_2_models_global_best_estim_1984-2013.png"), width=12,height=5,units="in",res=200)
mycols= colorRampPalette(brewer.pal(2,"Paired"))(2)
barplot2(burd_tot/1000, beside=T, plot.ci=T, ci.l = burd_tot_inf/1000, ci.u=burd_tot_sup/1000, col=mycols,
         ylab = "Severe cases (thousand)", ylim= c(0,850), las=2, main = "YF Burden, 1984-2013", cex.main=1.4, space = c(0,0.65))
legend("topright", legend=c( "Static model", "Dynamic model"), fill = mycols)
dev.off()



### by 5y band
ylist = list(rep(1:6, each=5))
burd_tot_foi_5y =aggregate(x=burd_tot_foi, by=ylist, FUN=sum)$x
burd_tot_foi_inf_5y=aggregate(x=burd_tot_foi_inf, by=ylist, FUN=sum)$x
burd_tot_foi_sup_5y=aggregate(x=burd_tot_foi_sup, by=ylist, FUN=sum)$x

burd_tot_R0_5y =aggregate(x=burd_tot_R0, by=ylist, FUN=sum)$x
burd_tot_R0_inf_5y=aggregate(x=burd_tot_R0_inf, by=ylist, FUN=sum)$x
burd_tot_R0_sup_5y=aggregate(x=burd_tot_R0_sup, by=ylist, FUN=sum)$x

burd_tot_5y = rbind(burd_tot_foi_5y, burd_tot_R0_5y)
burd_tot_inf_5y = rbind(burd_tot_foi_inf_5y, burd_tot_R0_inf_5y)
burd_tot_sup_5y = rbind(burd_tot_foi_sup_5y, burd_tot_R0_sup_5y)
colnames(burd_tot_5y) = colnames( burd_tot_inf_5y) = colnames( burd_tot_sup_5y) = c("1984-1988", "1989-1993", "1994-1998", "1999-2003", "2004-2008", "2009-2013")
#vec_years = c("1983-1987", "1988-1992", "1993-1997", "1998-2002", "2003-2007", "2008-2013")

png(paste0(outdir,"compare_burden_2_models_global_best_estim_1984-2013_by5y.png"), width=8,height=7,units="in",res=200)
mycols= colorRampPalette(brewer.pal(2,"Paired"))(2)
b = barplot2(burd_tot_5y/1000, beside=T, plot.ci=T, ci.l = burd_tot_inf_5y/1000, ci.u=burd_tot_sup_5y/1000, col=mycols,angle=45,
         ylab = "Severe cases (thousand)", ylim= c(0,4000), las=2, main = "YF Burden, 1984-2013", cex.main=1.4, space = c(0,0.65))
legend("topright", legend=c( "Static model", "Dynamic model"), fill = mycols)
dev.off()

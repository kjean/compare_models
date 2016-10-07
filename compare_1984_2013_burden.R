library(maptools)
library(sp) 
library(shapefiles)
library(fields)
library(Hmisc)
library(fields)

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




# import burden FOI model

foi_dir = paste(homedir, "script_MCMC/post_GAVI_runs_FOI_model/best_estimates/", sep="")
tab = read.csv(paste0(foi_dir,"cases_by_year_adm0_nb_runs=1000.csv"), h=T)
dim(tab)
colnames(tab)
head(tab[,1:4])
tab = tab[,-2]
colnames(tab)
tab[,2:ncol(tab)]=prop.death.all.cases*tab[,2:ncol(tab)] # apply CFR


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


# import 1984-2013 burden R0 model
setwd(paste(homedir, "script_MCMC_herd_immunity/post_GAVI_round/Thu_Mar_10_20164compartiments_finalMCMC_nb_runs=1000/best_estimates", sep=""))
tab = read.csv("cases_by_year_adm0_nb_runs=1000.csv", h=T)
dim(tab)
colnames(tab)
tab= tab[,-2]
tab[,2:ncol(tab)]=prop.death.all.cases*tab[,2:ncol(tab)] # apply CFR

burd_mean_R0 = burd_inf_R0 = burd_sup_R0 = NULL
for(adm in c34){
  print(adm)
  temp = tab[tab$adm0 == adm,]
  mean_tmp = apply(temp[,-1], 2, mean)
  low_tmp = apply(temp[,-1], 2, quantile, probs=0.025)
  sup_tmp = apply(temp[,-1], 2,  quantile, probs=0.975)
  burd_mean_R0 = rbind(burd_mean_R0, mean_tmp)
  burd_inf_R0 = rbind(burd_inf_R0, low_tmp)
  burd_sup_R0 = rbind(burd_sup_R0, sup_tmp)
}


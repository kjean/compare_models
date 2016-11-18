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
commondir = paste(homedir,"re-fit_models_2016/", sep="")
currdir = paste(homedir,"re-fit_models_2016/", "compare_models/", sep="")
outdir = paste0(currdir, "outputs_figures/")
outbreakdir = paste0(homedir, "Outbreak and vaccination/")

setwd(currdir)



scenar_type_vec = c("no_vaccine",  "response", "prev_campaigns_NoGAVI", "prev_campaigns", "routine_NoGAVI","routine" )

# scenar_type = "no_vaccine"
#scenar_type = "response"
#scenar_type = "prev_campaigns_NoGAVI"
#scenar_type = "prev_campaigns"
#scenar_type = "routine_NoGAVI"
#scenar_type = "routine"

computer<-2
if (computer==0) { homedir<-"Y:/Kevin/"
} else if (computer==1) { homedir<-"C:" 
} else if (computer==2) { homedir<-"//fi--didenas1/yf/Kevin/"}


shpdir = paste(homedir,"shapefiles/gadm2/",sep="")
currdir = paste(homedir,"GAVI 2015 estimates/", sep="")


c34 = c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ERI", "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GNQ", "KEN", "LBR", "MLI", "MRT", "NER", "NGA", "RWA", "SDN", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZMB")
country34 = c("Angola","Burundi","Benin","Burkina Faso","Central African Rep.","C?te d'Ivoire","Cameroon","Dem. Rep. of Congo","Rep. of Congo","Eritrea","Ethiopia","Gabon","Ghana","Guinea","Gambia","Guinea-Bissau","Equatorial Guinea","Kenya","Liberia","Mali","Mauritania","Niger","Nigeria","Rwanda","Sudan","Senegal","Sierra Leone","Somalia", "South Sudan", "Chad","Togo","Tanzania","Uganda", "Zambia") 


set.seed(101)
prop.death.all.cases = rlnorm(1000, meanlog= -3.00, sdlog=0.442)
prop.death.all.cases = rep(prop.death.all.cases, each = 34)

########################################
#### VACCINATION IMPACT
########################################
#### 5 and 95percentil of vaccination impact 2001-2030




###########  ###########  ###########  ###########  ###########  
###########  function
###########  ###########  ###########  ###########  ###########  

calculate_vaccine_impact = function(scenar, global_country="country", model){
  
  scenar_type_vec = c("no_vaccine",  "response", "prev_campaigns_NoGAVI", "prev_campaigns", "routine_NoGAVI","routine" )
  if( !(scenar %in% scenar_type_vec) ){
    stop("Wrong scenar type; choose in <no_vaccine>,  <response>, <prev_campaigns_NoGAVI>, <prev_campaigns>, <routine_NoGAVI>,<routine> ")
  }
  
  if( !(global_country %in% c("country", "global") ) ){
    stop("Wrong global_country argument: choose <global> or <country> ")
  }
  
  if( !(model %in% c("FOI", "R0") ) ){stop("Wrong model specification: model has to be <FOI> or <R0>")}
  
  set.seed(101)
  prop.death.all.cases = rlnorm(1000, meanlog= -3.00, sdlog=0.442)
  prop.death.all.cases = rep(prop.death.all.cases, each = 34)
  
  if(scenar == "no_vaccine"){year_cut1 = 2000 
  } else if(scenar == "response"){
    year_cut1 = 2017
  } else if(scenar == "prev_campaigns_NoGAVI" | scenar =="prev_campaigns"){
    year_cut1 = 2023
  } else if(scenar == "routine_NoGAVI" | scenar =="routine"){
    year_cut1 = 2030
  }
  
  if(model=="R0"){
    tab = read.csv(paste(currdir, "R0_model/",scenar, "/all_cases_", scenar, "_", year_cut1, "_nb_runs=1000_by_year_R0", ".csv", sep=""), h=T)
  } else if(model=="FOI"){
    tab = read.csv(paste(currdir, "FOI_model/",scenar, "/all_cases_", scenar, "_", year_cut1, "_nb_runs=1000_by_year_", ".csv", sep=""), h=T)
  }
  
  tab[,2:ncol(tab)]=prop.death.all.cases*tab[,2:ncol(tab)]
  nbruns = rep(1:1000, each=34)
  
  if(global_country == "global"){
    burden = NULL
    for(i in 1:1000){
      sub_tab = tab[nbruns ==i ,2:31]
      burden_tmp = sum(sub_tab)
      burden = c(burden, burden_tmp)
      
    }
    
    burden_mean = mean(burden)
    burden_median = median(burden)
    
    burden_95_low = quantile(burden, prob = 0.025)
    burden_95_sup = quantile(burden, prob = 0.975)
    burden_90_low = quantile(burden, prob = 0.05)
    burden_90_sup = quantile(burden, prob = 0.095)
    
  } else if(global_country == "country"){
    burden = NULL
    for(i in 1:1000){
      sub_tab = tab[nbruns ==i ,2:31]
      burden_tmp = rowSums(sub_tab)
      burden = cbind(burden, burden_tmp)
    }
    burden_mean = apply(burden, 1, mean)
    burden_median = apply(burden, 1, median)
    
    burden_95_low = apply(burden, 1, quantile, prob = 0.025)
    burden_95_sup = apply(burden, 1, quantile, prob = 0.975)
    burden_90_low = apply(burden, 1, quantile, prob = 0.05)
    burden_90_sup = apply(burden, 1, quantile, prob = 0.95)
    
    names(burden_mean) = names(burden_median) = names(burden_95_low) = names(burden_95_sup) = names(burden_90_low) = names(burden_90_sup) = c34
    
  }
  
  burden_list = list(burden_mean = burden_mean, burden_median=burden_median, burden_95_low=burden_95_low, burden_95_sup=burden_95_sup,
                     burden_90_low=burden_90_low, burden_90_sup=burden_90_sup)  
  
  return(burden_list)
}



########
# same function but can specify years of calculation
calculate_vaccine_impact_specify_years = function(scenar, global_country="country", model, years){
  
  scenar_type_vec = c("no_vaccine",  "response", "prev_campaigns_NoGAVI", "prev_campaigns", "routine_NoGAVI","routine" )
  if( !(scenar %in% scenar_type_vec) ){
    stop("Wrong scenar type; choose in <no_vaccine>,  <response>, <prev_campaigns_NoGAVI>, <prev_campaigns>, <routine_NoGAVI>,<routine> ")
  }
  
  if( !(global_country %in% c("country", "global") ) ){
    stop("Wrong global_country argument: choose <global> or <country> ")
  }
  
  if( !(model %in% c("FOI", "R0") ) ){stop("Wrong model specification: model has to be <FOI> or <R0>")}
  
  set.seed(101)
  prop.death.all.cases = rlnorm(1000, meanlog= -3.00, sdlog=0.442)
  prop.death.all.cases = rep(prop.death.all.cases, each = 34)
  
  if(scenar == "no_vaccine"){year_cut1 = 2000 
  } else if(scenar == "response"){
    year_cut1 = 2017
  } else if(scenar == "prev_campaigns_NoGAVI" | scenar =="prev_campaigns"){
    year_cut1 = 2023
  } else if(scenar == "routine_NoGAVI" | scenar =="routine"){
    year_cut1 = 2030
  }
  
  if(model=="R0"){
    tab = read.csv(paste(currdir, "R0_model/",scenar, "/all_cases_", scenar, "_", year_cut1, "_nb_runs=1000_by_year_R0", ".csv", sep=""), h=T)
  } else if(model=="FOI"){
    tab = read.csv(paste(currdir, "FOI_model/",scenar, "/all_cases_", scenar, "_", year_cut1, "_nb_runs=1000_by_year_", ".csv", sep=""), h=T)
  }
  colnames(tab) = c("adm0", 2001:2100)
  mm.years = match(years, colnames(tab))
  
  tab=prop.death.all.cases*tab[,mm.years]
  nbruns = rep(1:1000, each=34)
  
  if(global_country == "global"){
    burden = NULL
    for(i in 1:1000){
      sub_tab = tab[nbruns ==i ,1:ncol(tab)]
      burden_tmp = sum(sub_tab)
      burden = c(burden, burden_tmp)
      
    }
    
    burden_mean = mean(burden)
    burden_median = median(burden)
    
    burden_95_low = quantile(burden, prob = 0.025)
    burden_95_sup = quantile(burden, prob = 0.975)
    burden_90_low = quantile(burden, prob = 0.05)
    burden_90_sup = quantile(burden, prob = 0.095)
    
  } else if(global_country == "country"){
    burden = NULL
    for(i in 1:1000){
      sub_tab = tab[nbruns ==i ,1:ncol(tab)]
      burden_tmp = rowSums(sub_tab)
      burden = cbind(burden, burden_tmp)
    }
    burden_mean = apply(burden, 1, mean)
    burden_median = apply(burden, 1, median)
    
    burden_95_low = apply(burden, 1, quantile, prob = 0.025)
    burden_95_sup = apply(burden, 1, quantile, prob = 0.975)
    burden_90_low = apply(burden, 1, quantile, prob = 0.05)
    burden_90_sup = apply(burden, 1, quantile, prob = 0.95)
    
    names(burden_mean) = names(burden_median) = names(burden_95_low) = names(burden_95_sup) = names(burden_90_low) = names(burden_90_sup) = c34
    
  }
  
  burden_list = list(burden_mean = burden_mean, burden_median=burden_median, burden_95_low=burden_95_low, burden_95_sup=burden_95_sup,
                     burden_90_low=burden_90_low, burden_90_sup=burden_90_sup)  
  
  return(burden_list)
}





##### calculing and plotting
scenar_type_vec = c("no_vaccine",  "response", "prev_campaigns_NoGAVI", "prev_campaigns", "routine_NoGAVI","routine" )



#### R0 model
burden_no_vaccine_R0 = calculate_vaccine_impact("no_vaccine", "global", "R0")
burden_response_R0 = calculate_vaccine_impact("response", "global", "R0")
burden_prev_campaigns_R0 = calculate_vaccine_impact("prev_campaigns", "global", "R0")
burden_routine_R0 = calculate_vaccine_impact("routine", "global", "R0")



impact_response_mean_R0 = burden_no_vaccine_R0$burden_mean - burden_response_R0$burden_mean
impact_response_95_low_R0 = burden_no_vaccine_R0$burden_95_low - burden_response_R0$burden_95_low
impact_response_95_sup_R0 = burden_no_vaccine_R0$burden_95_sup - burden_response_R0$burden_95_sup

impact_prev_campaigns_mean_R0 = burden_response_R0$burden_mean - burden_prev_campaigns_R0$burden_mean
impact_prev_campaigns_95_low_R0 = burden_response_R0$burden_95_low - burden_prev_campaigns_R0$burden_95_low
impact_prev_campaigns_95_sup_R0 = burden_response_R0$burden_95_sup - burden_prev_campaigns_R0$burden_95_sup

impact_routine_mean_R0 = burden_prev_campaigns_R0$burden_mean - burden_routine_R0$burden_mean
impact_routine_95_low_R0 = burden_prev_campaigns_R0$burden_95_low - burden_routine_R0$burden_95_low
impact_routine_95_sup_R0 = burden_prev_campaigns_R0$burden_95_sup - burden_routine_R0$burden_95_sup

burden_mean_R0 = c(burden_no_vaccine_R0$burden_mean, burden_response_R0$burden_mean,
                   burden_prev_campaigns_R0$burden_mean, burden_routine_R0$burden_mean)
burden_95_low_R0 = c(burden_no_vaccine_R0$burden_95_low, burden_response_R0$burden_95_low,
                     burden_prev_campaigns_R0$burden_95_low, burden_routine_R0$burden_95_low)
burden_95_sup_R0 = c(burden_no_vaccine_R0$burden_95_sup, burden_response_R0$burden_95_sup,
                     burden_prev_campaigns_R0$burden_95_sup, burden_routine_R0$burden_95_sup)

impact_R0 = c(impact_response_mean_R0, impact_prev_campaigns_mean_R0, impact_routine_mean_R0)
impact_R0_low = c(impact_response_95_low_R0, impact_prev_campaigns_95_low_R0, impact_routine_95_low_R0)
impact_R0_sup = c(impact_response_95_sup_R0, impact_prev_campaigns_95_sup_R0, impact_routine_95_sup_R0)


#### FOI model
burden_no_vaccine_FOI = calculate_vaccine_impact("no_vaccine", "global", "FOI")
burden_response_FOI = calculate_vaccine_impact("response", "global", "FOI")
burden_prev_campaigns_FOI = calculate_vaccine_impact("prev_campaigns", "global", "FOI")
burden_routine_FOI = calculate_vaccine_impact("routine", "global", "FOI")



impact_response_mean_FOI = burden_no_vaccine_FOI$burden_mean - burden_response_FOI$burden_mean
impact_response_95_low_FOI = burden_no_vaccine_FOI$burden_95_low - burden_response_FOI$burden_95_low
impact_response_95_sup_FOI = burden_no_vaccine_FOI$burden_95_sup - burden_response_FOI$burden_95_sup

impact_prev_campaigns_mean_FOI = burden_response_FOI$burden_mean - burden_prev_campaigns_FOI$burden_mean
impact_prev_campaigns_95_low_FOI = burden_response_FOI$burden_95_low - burden_prev_campaigns_FOI$burden_95_low
impact_prev_campaigns_95_sup_FOI = burden_response_FOI$burden_95_sup - burden_prev_campaigns_FOI$burden_95_sup

impact_routine_mean_FOI = burden_prev_campaigns_FOI$burden_mean - burden_routine_FOI$burden_mean
impact_routine_95_low_FOI = burden_prev_campaigns_FOI$burden_95_low - burden_routine_FOI$burden_95_low
impact_routine_95_sup_FOI = burden_prev_campaigns_FOI$burden_95_sup - burden_routine_FOI$burden_95_sup

burden_mean_FOI = c(burden_no_vaccine_FOI$burden_mean, burden_response_FOI$burden_mean,
                   burden_prev_campaigns_FOI$burden_mean, burden_routine_FOI$burden_mean)
burden_95_low_FOI = c(burden_no_vaccine_FOI$burden_95_low, burden_response_FOI$burden_95_low,
                     burden_prev_campaigns_FOI$burden_95_low, burden_routine_FOI$burden_95_low)
burden_95_sup_FOI = c(burden_no_vaccine_FOI$burden_95_sup, burden_response_FOI$burden_95_sup,
                     burden_prev_campaigns_FOI$burden_95_sup, burden_routine_FOI$burden_95_sup)

impact_FOI = c(impact_response_mean_FOI, impact_prev_campaigns_mean_FOI, impact_routine_mean_FOI)
impact_FOI_low = c(impact_response_95_low_FOI, impact_prev_campaigns_95_low_FOI, impact_routine_95_low_FOI)
impact_FOI_sup = c(impact_response_95_sup_FOI, impact_prev_campaigns_95_sup_FOI, impact_routine_95_sup_FOI)


burden_by_strat = cbind(burden_mean_FOI, burden_mean_R0)/1000000
burden_by_strat_low = cbind(burden_95_low_FOI, burden_95_low_R0)/1000000
burden_by_strat_sup = cbind(burden_95_sup_FOI, burden_95_sup_R0)/1000000
colnames(burden_by_strat) = c("Static model", "Dynamic model")



impact_by_strat = cbind(impact_FOI, impact_R0)/1000000
impact_by_strat_low = cbind(impact_FOI_low, impact_R0_low)/1000000
impact_by_strat_sup = cbind(impact_FOI_sup, impact_R0_sup)/1000000
colnames(impact_by_strat) = c("Static model", "Dynamic model")



response_doses = 30756758
prev_campaigns_doses = 186179404
routine_doses = 468846483

impact_by_strat_doses = impact_by_strat
impact_by_strat_doses[1,] = 1000000*1000*impact_by_strat[1,]/response_doses
impact_by_strat_doses[2,] = 1000000*1000*impact_by_strat[2,]/prev_campaigns_doses
impact_by_strat_doses[3,] = 1000000*1000*impact_by_strat[3,]/routine_doses

impact_by_strat_doses_low = impact_by_strat_low
impact_by_strat_doses_low[1,] = 1000000*1000*impact_by_strat_low[1,]/response_doses
impact_by_strat_doses_low[2,] = 1000000*1000*impact_by_strat_low[2,]/prev_campaigns_doses
impact_by_strat_doses_low[3,] = 1000000*1000*impact_by_strat_low[3,]/routine_doses

impact_by_strat_doses_sup = impact_by_strat_sup
impact_by_strat_doses_sup[1,] = 1000000*1000*impact_by_strat_sup[1,]/response_doses
impact_by_strat_doses_sup[2,] = 1000000*1000*impact_by_strat_sup[2,]/prev_campaigns_doses
impact_by_strat_doses_sup[3,] = 1000000*1000*impact_by_strat_sup[3,]/routine_doses



png("compare_burden_and_vaccine_impact_2_models_global.png", width=12,height=5,units="in",res=200)

par(mfrow = c(1,3))
#display.brewer.all()
mycols= colorRampPalette(brewer.pal(4,"Paired"))(nrow(burden_by_strat))
barplot2(burden_by_strat, beside=T, plot.ci=T, ci.l = burden_by_strat_low, ci.u=burden_by_strat_sup, col=mycols,
         ylab = "Deaths (millions)", main = "A. Cumulative 2001-2030 burden", cex.lab = 1.3, space = c(0,0.65),
         cex.main=1.4)
legend(1, 24, legend=c("No vaccine", "Reactive campaigns", "+ Preventive campaigns", "+ Routine"), fill = mycols)
#dev.off()


mycols= colorRampPalette(brewer.pal(4,"Paired"))(nrow(burden_by_strat))[2:4]
barplot2(impact_by_strat, beside=T, plot.ci=T, ci.l = impact_by_strat_low, ci.u=impact_by_strat_sup, col=mycols,
         ylab = "Deaths averted (millions)", main = "B. Vaccine impact (2001-2030)", cex.lab = 1.3, space = c(0,0.65),
         cex.main=1.4)
legend(1, 11, legend=c( "Reactive campaigns", "Preventive campaigns", "Routine"), fill = mycols)



mycols= colorRampPalette(brewer.pal(4,"Paired"))(nrow(burden_by_strat))[2:4]
barplot2(impact_by_strat_doses, beside=T, plot.ci=T, ci.l = impact_by_strat_doses_low, ci.u=impact_by_strat_doses_sup, col=mycols,
         ylab = "Deaths averted for 1,000 doses", main = "C. Vaccine impact per doses (2001-2030)", cex.lab = 1.3, space = c(0,0.65),
         cex.main=1.4)
legend(1, 33, legend=c( "Reactive campaigns", "Preventive campaigns", "Routine"), fill = mycols)

dev.off()


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### By country - total vaccine impact
burden_no_vaccine_R0 = calculate_vaccine_impact("no_vaccine", "country", "R0")
burden_response_R0 = calculate_vaccine_impact("response", "country", "R0")
burden_prev_campaigns_R0 = calculate_vaccine_impact("prev_campaigns", "country", "R0")
burden_routine_R0 = calculate_vaccine_impact("routine", "country", "R0")

impact_all_vac_mean_R0 = burden_no_vaccine_R0$burden_mean - burden_routine_R0$burden_mean
impact_all_vac_median_R0 = burden_no_vaccine_R0$burden_median - burden_routine_R0$burden_median
impact_all_vac_95_low_R0 = burden_no_vaccine_R0$burden_95_low - burden_routine_R0$burden_95_low
impact_all_vac_95_sup_R0 = burden_no_vaccine_R0$burden_95_sup - burden_routine_R0$burden_95_sup


burden_no_vaccine_FOI = calculate_vaccine_impact("no_vaccine", "country", "FOI")
burden_response_FOI = calculate_vaccine_impact("response", "country", "FOI")
burden_prev_campaigns_FOI = calculate_vaccine_impact("prev_campaigns", "country", "FOI")
burden_routine_FOI = calculate_vaccine_impact("routine", "country", "FOI")

impact_all_vac_mean_FOI = burden_no_vaccine_FOI$burden_mean - burden_routine_FOI$burden_mean
impact_all_vac_median_FOI = burden_no_vaccine_FOI$burden_median - burden_routine_FOI$burden_median
impact_all_vac_95_low_FOI = burden_no_vaccine_FOI$burden_95_low - burden_routine_FOI$burden_95_low
impact_all_vac_95_sup_FOI = burden_no_vaccine_FOI$burden_95_sup - burden_routine_FOI$burden_95_sup

impact_all_vac_mean = rbind(impact_all_vac_mean_FOI, impact_all_vac_mean_R0)
impact_all_vac_median = rbind(impact_all_vac_median_FOI, impact_all_vac_median_R0)
impact_all_vac_95_low = rbind(impact_all_vac_95_low_FOI, impact_all_vac_95_low_FOI)
impact_all_vac_95_sup = rbind(impact_all_vac_95_sup_FOI, impact_all_vac_95_sup_R0 )

barplot2(impact_all_vac_mean, beside=T, col=c("grey50", "grey80"))


# try to group by country
colnames(impact_all_vac_mean )

barplot2(impact_all_vac_mean[,-which(colnames(impact_all_vac_mean)=="NGA")], beside=T, col=c("grey50", "grey80"))

country_to_retrieve = c("BDI", "ERI", "GAB", "GNQ", "RWA", "SOM", "TZA", "ZMB")

c23 = c("AGO","BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GMB", "GNB", "KEN", "LBR", "MLI", "MRT", "NER",  
        "SDN", "SEN", "SLE", "SSD", "TCD", "TGO", "UGA")
country23 = country34[c34 %in% c23]


impact_all_vac_mean_23 = impact_all_vac_mean[, which(colnames(impact_all_vac_mean) %in% c23)]
impact_all_vac_95_low_23 = impact_all_vac_95_low[, which(colnames(impact_all_vac_95_low) %in% c23)]
impact_all_vac_95_sup_23 = impact_all_vac_95_sup[, which(colnames(impact_all_vac_95_sup) %in% c23)]
colnames(impact_all_vac_mean_23) = colnames(impact_all_vac_95_low_23) =  colnames(impact_all_vac_95_sup_23) = country23


impact_all_vac_mean_NGA = matrix(impact_all_vac_mean[, which(colnames(impact_all_vac_mean)=="NGA")])
colnames(impact_all_vac_mean_NGA) = "NGA"
impact_all_vac_95_low_NGA = impact_all_vac_95_low[, which(colnames(impact_all_vac_95_low) =="NGA")]
impact_all_vac_95_sup_NGA = impact_all_vac_95_sup[, which(colnames(impact_all_vac_95_sup) =="NGA")]
colnames(impact_all_vac_mean_NGA) =  "Nigeria"

setwd(currdir)
png("compare_burden_2_models_countries.png", width=12,height=6,units="in",res=100)
layout(matrix(c(rep(1,9),2), nrow=1 ))
par(oma= c(10,4,4,4))
par(mgp=c(5,1,0))
par(mar=c(4,7,4,0))
barplot2(impact_all_vac_mean_23, beside=T, plot.ci=T, ci.l = impact_all_vac_95_low_23, ci.u=impact_all_vac_95_sup_23, col=c("grey50", "grey80"),
         ylab = "Number of deaths averted (2001-2030)", cex.lab = 1.5, cex.names = 1.4, ci.width = 0.3, ylim=c(0, 1300000), 
         las=2, space = c(0,0.65))
par(mar=c(4,2,4,4))
barplot2(impact_all_vac_mean_NGA, beside=T, plot.ci=T, ci.l = impact_all_vac_95_low_NGA, ci.u=impact_all_vac_95_sup_NGA, col=c("grey50", "grey80"),
         ylab = "Nb of deaths averted", cex.names = 1.4, ci.width = 0.2, las=2)
dev.off()

# 
# 
# barplot2(burd_med_log, beside=T, plot.ci=T, ci.l = burd_low_log, ci.u=burd_sup_log, col=c("grey50", "grey80"), 
#          names.arg=c("Model without herd effect","Model with herd effect"), legend.text=c("2005", "2013") ,
#          ci.lwd = 1.5, ci.width = 0.1, ylim=c(3, max(burd_sup_log)), xpd=F,
#          axes=F, ylab = "Number of deaths",
#          main = "Overall burden estimate comparison, log scale")
# 
# 
# 




############################
# Compare burden on specific years
###########################

burden_2006_2013_routine_R0 = calculate_vaccine_impact_specify_years("routine", global_country="global", "R0", years=2006:2013)
burden_2006_2013_routine_foi = calculate_vaccine_impact_specify_years("routine", global_country="global", "FOI", years=2006:2013)

burden_2006_2013_no_vaccine_R0 = calculate_vaccine_impact_specify_years("no_vaccine", global_country="global", "R0", years=2006:2013)
burden_2006_2013_no_vaccine_foi = calculate_vaccine_impact_specify_years("no_vaccine", global_country="global", "FOI", years=2006:2013)

impact_vaccin_foi_mean = burden_2006_2013_no_vaccine_foi$burden_mean - burden_2006_2013_routine_foi$burden_mean
impact_vaccin_foi_inf = burden_2006_2013_no_vaccine_foi$burden_95_low - burden_2006_2013_routine_foi$burden_95_low
impact_vaccin_foi_sup = burden_2006_2013_no_vaccine_foi$burden_95_sup - burden_2006_2013_routine_foi$burden_95_sup

impact_vaccin_R0_mean = burden_2006_2013_no_vaccine_R0$burden_mean - burden_2006_2013_routine_R0$burden_mean
impact_vaccin_R0_inf = burden_2006_2013_no_vaccine_R0$burden_95_low - burden_2006_2013_routine_R0$burden_95_low
impact_vaccin_R0_sup = burden_2006_2013_no_vaccine_R0$burden_95_sup - burden_2006_2013_routine_R0$burden_95_sup



## number of doses from the excel sheets sent to GAVI - need to consider all doses from 2001 to 2013 as here are the difference 
# could be decided otherwise though...
response_doses = 30756758
prev_campaigns_doses = 70251779.8
routine_doses = 106119990

nb_doses = response_doses + prev_campaigns_doses + routine_doses

1000*impact_vaccin_foi_mean /nb_doses
1000*impact_vaccin_foi_inf /nb_doses
1000*impact_vaccin_foi_sup /nb_doses

1000*impact_vaccin_R0_mean /nb_doses
1000*impact_vaccin_R0_inf /nb_doses
1000*impact_vaccin_R0_sup /nb_doses


# impact all vaccination
impact = burden_no_vaccine - burden_routine
hist(impact)
median(impact)
mean(impact)
quantile(impact, probs= c(0.05, 0.95))


# impact GAVI

# impact routine GAVI
impact_routine_GAVI = burden_routine_NoGAVI - burden_routine
mean(impact_routine_GAVI)

# impact campaigns GAVI
impact_camp_GAVI = burden_prev_campaigns_NoGAVI - burden_prev_campaigns
mean(impact_camp_GAVI)
  
# impact response
impact_response = burden_no_vaccine - burden_response
mean(impact_response)


impact_GAVI = impact_routine_GAVI + impact_camp_GAVI + impact_response
mean(impact_GAVI)
quantile(impact_GAVI, probs= c(0.05, 0.95))



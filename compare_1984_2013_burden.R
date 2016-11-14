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



####### ####### ####### ####### ####### 
####### load shapefiles


# shp0 = readShapePoly(paste(shpdir, "Africa_adm0.shp",sep=""))
# shp1 = readShapePoly(paste(shpdir, "Africa_adm1.shp",sep=""))
# 
# shp1$adm0_adm1<-paste(shp1$ISO, shp1$ID_1, sep="_")
# shp1 = shp1[order(shp1$adm0_adm1),]

c34 = c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ERI", "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GNQ", "KEN", "LBR", "MLI", "MRT", "NER", "NGA", "RWA", "SDN", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZMB")
country34 = c("Angola","Burundi","Benin","Burkina Faso","Central African Republic","Côte d'Ivoire","Cameroon","Democratic Republic of the Congo","Republic of Congo","Eritrea","Ethiopia","Gabon","Ghana","Guinea","Gambia","Guinea-Bissau","Equatorial Guinea","Kenya","Liberia","Mali","Mauritania","Niger","Nigeria","Rwanda","Sudan","Senegal","Sierra Leone","Somalia", "South Sudan", "Chad","Togo","Tanzania","Uganda", "Zambia") 


################ outbreak data
#### read outbreak data
outbreak = read.csv(paste0(outbreakdir, "outbreaks_1960s-2014_KJ.csv"), h=T, stringsAsFactors = F)
str(outbreak)
outbreak$year = as.numeric(outbreak$year)
table(outbreak$year, exclude = NULL)
outbreak = outbreak[outbreak$year >1983 & outbreak$year<2014,]
length(table(outbreak$year))
hist(outbreak$year)








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

foi_dir = paste0(commondir, "YF_FOI_burden_model/")
tab = read.csv(paste0(foi_dir,"all_cases_FOI_best_estimate_nb_runs=1000_by_year_.csv"), h=T)
dim(tab)
colnames(tab)
head(tab[,1:4])
#tab = tab[,-2]
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
r0_dir =  paste0(commondir, "YF_R0_burden_model/")
tab = read.csv(paste0(r0_dir,"all_cases_R0_best_estimate_nb_runs=1000_by_year_.csv"), h=T)
dim(tab)
colnames(tab)
#tab= tab[,-2]
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
#colnames(burd_tot_5y) = colnames( burd_tot_inf_5y) = colnames( burd_tot_sup_5y) = c("1984-1988", "1989-1993", "1994-1998", "1999-2003", "2004-2008", "2009-2013")
vec_years = c("1984-1988", "1989-1993", "1994-1998", "1999-2003", "2004-2008", "2009-2013")


# outbreak data by 5y band
outbreak$year_5 = cut(outbreak$year, breaks = c(1983,1988, 1993, 1998, 2003, 2008, 2013))
table(outbreak$year_5)
tab5 = t(matrix(table(outbreak$year_5)))


png(paste0(outdir,"compare_burden_2_models_1984-2013_by5y_w_outbreaks.png"), width=9,height=5,units="in",res=200)
par(oma=c(0,4,0,4))
mycols= colorRampPalette(brewer.pal(2,"Paired"))(2)
b = barplot2(burd_tot_5y/1000, beside=T, plot.ci=T, ci.l = burd_tot_inf_5y/1000, ci.u=burd_tot_sup_5y/1000, col=mycols,angle=45,
         ylab = "", ylim= c(0,4000), las=2, main = "YF Burden, 1984-2013", cex.main=1.4, space = c(0,0.65))
text(apply(b, 2, mean)+0.5, par("usr")[3]-100, srt = 45, adj= 1, xpd = TRUE,
     labels = vec_years, cex=1)
tab5_rescalled = tab5*3500/118
points(apply(b, 2, mean), tab5_rescalled, type="l", col="tomato", lwd=2 )
points(apply(b, 2, mean), tab5_rescalled, pch=19, col="tomato", lwd=2)
ticks_ax_2 = c(0,25,50,75,100,125)
ticks_ax_2_scaled = ticks_ax_2*3500/118
axis(side=4, at=ticks_ax_2_scaled, labels=FALSE, col="tomato", lwd=2)
text(x =par("usr")[2]+0.7, y = ticks_ax_2_scaled,  labels = ticks_ax_2,col="tomato", xpd=T)
legend(par("usr")[1]+1, par("usr")[4]-100, legend=c( "Static model", "Dynamic model"), fill = mycols)
mtext(side=2, "Estimated nb of severe cases", outer=T,cex=1.3, font=2)
mtext(side=4, "Number of outbreaks reported", outer=T, line=0.7, cex=1.3, font=2, col="tomato")
dev.off()



### by 3y band
ylist = list(rep(1:10, each=3))
burd_tot_foi_3y =aggregate(x=burd_tot_foi, by=ylist, FUN=sum)$x
burd_tot_foi_inf_3y=aggregate(x=burd_tot_foi_inf, by=ylist, FUN=sum)$x
burd_tot_foi_sup_3y=aggregate(x=burd_tot_foi_sup, by=ylist, FUN=sum)$x

burd_tot_R0_3y =aggregate(x=burd_tot_R0, by=ylist, FUN=sum)$x
burd_tot_R0_inf_3y=aggregate(x=burd_tot_R0_inf, by=ylist, FUN=sum)$x
burd_tot_R0_sup_3y=aggregate(x=burd_tot_R0_sup, by=ylist, FUN=sum)$x

burd_tot_3y = rbind(burd_tot_foi_3y, burd_tot_R0_3y)
burd_tot_inf_3y = rbind(burd_tot_foi_inf_3y, burd_tot_R0_inf_3y)
burd_tot_sup_3y = rbind(burd_tot_foi_sup_3y, burd_tot_R0_sup_3y)
colnames(burd_tot_3y) = colnames( burd_tot_inf_3y) = colnames( burd_tot_sup_3y) = c("1984-1986", "1987-1989", "1990-1992", "1993-1995", "1996-1998", "1999-2001","2002-2004", "2005-2007", "2008-2010", "2011-2013")
#vec_years = c("1983-1987", "1988-1992", "1993-1997", "1998-2002", "2003-2007", "2008-2013")

png(paste0(outdir,"compare_burden_2_models_global_best_estim_1984-2013_by3y.png"), width=8,height=7,units="in",res=200)
mycols= colorRampPalette(brewer.pal(2,"Paired"))(2)
b = barplot2(burd_tot_3y/1000, beside=T, plot.ci=T, ci.l = burd_tot_inf_3y/1000, ci.u=burd_tot_sup_3y/1000, col=mycols,angle=45,
             ylab = "Severe cases (thousand)", ylim= c(0,4000), las=2, main = "YF Burden, 1984-2013", cex.main=1.4, space = c(0,0.65))
legend("topright", legend=c( "Static model", "Dynamic model"), fill = mycols)
dev.off()
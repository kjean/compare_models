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



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
outbreakdir = paste0(homedir, "Outbreak and vaccination/")

setwd(currdir)


#### read outbreak data
outbreak = read.csv(paste0(outbreakdir, "outbreaks_1960s-2014_KJ.csv"), h=T, stringsAsFactors = F)
str(outbreak)
outbreak$year = as.numeric(outbreak$year)
table(outbreak$year, exclude = NULL)
outbreak = outbreak[outbreak$year >1983 & outbreak$year<2014,]
length(table(outbreak$year))
hist(outbreak$year)

# by 5y band
outbreak$year_3 = cut(outbreak$year, breaks = c(1983,1986, 1989, 1992, 1995, 1998, 2001, 2004, 2007, 2010, 2013))
table(outbreak$year_3)

tab5 = t(matrix(table(outbreak$year_5)))
colnames( tab5) = c("1984-1988", "1989-1993", "1994-1998", "1999-2003", "2004-2008", "2009-2013")
mycol= colorRampPalette(brewer.pal(4,"Paired"))(4)[4]
barplot(tab5, col = mycol,  ylab = "Number of outbreaks reported", width=0.2, space = 0.2)


# by 3y band
outbreak$year_3 = cut(outbreak$year, breaks = c(1983,1988, 1993, 1998, 2003, 2008, 2013))
table(outbreak$year_3)

tab3 = t(matrix(table(outbreak$year_3)))
colnames( tab3) = c("1984-1986", "1987-1989", "1990-1992", "1993-1995", "1996-1998", "1999-2001","2002-2004", "2005-2007", "2008-2010", "2011-2013")
mycol= colorRampPalette(brewer.pal(4,"Paired"))(4)[4]
barplot(tab3, col = mycol,  ylab = "Number of outbreaks reported", width=0.2, space = 0.2)
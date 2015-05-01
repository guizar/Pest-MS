#/////////
# IPM SUMMARIES
# evaluate all code and look for the following objects:

# SUM_IPM: IPM Summary across all temp anomalies and phis
# SUM_IPM_REGIONAL: IPM Summary across all temp anomalies, phis and regions
# SUM_IPM_TROPICS: IPM Summary across all temp anomalies and phis within the tropics (absolute LAT < 23)
#////////

rm(list = ls())

# Load data and initial variables ------

wd = "~/R/Pest-MS/"
wdpng = "~/R/Pest-MS/png"
wdtables = "~/R/Pest-MS/tables"
wddata = "~/R/Pest-MS/data/"
wdrdata = "~/R/Pest-MS/RData/"
wdfun = "~/R/Pest-MS/fun"

# Load data
load(file.path(wdrdata,"ALL_2c.RData"))
load(file.path(wdrdata,"ALL_4c.RData"))

# clean space
# source(file.path(wdfun,"clean_space.r")) 

# libs
library(ggplot2)
library(plyr)
library(dplyr)

# functions
source(file.path(wdfun,"cols_gg.r")) 
source(file.path(wdfun,"multiplots_gg.r"))
source(file.path(wdfun,"custom_cut.r"))

# --- 1) REMOVE IPM OUTLIERS FROM MASTER TABLES ----
DAT_2c = ddply(DAT_2c,c("fact","phi","crop"), function(d){
  limits = median(d$value) + 3*c(-1, 1)*mad(d$value)
  d[(d$value - limits[1])*(limits[2] - d$value) > 0,]
})

DAT_4c = ddply(DAT_4c,c("fact","phi","crop"), function(d){
  limits = median(d$value) + 3*c(-1, 1)*mad(d$value)
  d[(d$value - limits[1])*(limits[2] - d$value) > 0,]
})


#--- 2) SUMMARY BY TEMP ANOM (2c and 4c) ----
SUM_2c = ddply(DAT_2c,c("fact","phi"),summarise,
           mean= mean(value,na.rm = T),
           median= median(value,na.rm = T),
           min= min(value,na.rm = T),
           max= max(value,na.rm = T))
SUM_2c$anomaly = "2C"

SUM_4c = ddply(DAT_4c,c("fact","phi"),summarise,
           mean= mean(value,na.rm = T),
           median= median(value,na.rm = T),
           min= min(value,na.rm = T),
           max= max(value,na.rm = T))
SUM_4c$anomaly = "4C"

# SUMMARY IPM
SUM = rbind(SUM_2c,SUM_4c)
SUM_IPM = dplyr::filter(SUM,fact=="IPM")


#--- 3) REGIONAL IPM SUMMARY (2c and 4c) ----
SUM_2c = ddply(DAT_2c,c("fact","phi","region"),summarise,
              mean= mean(value,na.rm = T),
              median= median(value,na.rm = T),
              min= min(value,na.rm = T),
              max= max(value,na.rm = T))
SUM_2c$anomaly = "2C"

SUM_4c = ddply(DAT_4c,c("fact","phi","region"),summarise,
              mean= mean(value,na.rm = T),
              median= median(value,na.rm = T),
              min= min(value,na.rm = T),
              max= max(value,na.rm = T))
SUM_4c$anomaly = "4C"

# SUMMARY IPM
SUM = rbind(SUM_2c,SUM_4c)
SUM_IPM_REGIONAL = dplyr::filter(SUM,fact=="IPM")


rm(SUM_2c)
rm(SUM_4c)
rm(SUM)


#--- 4) TROPICAL IPM SUMMARY (2c and 4c) ----

# --- 2c --
fg = dplyr::select(ALL_2c, one_of("LAT",grep(pattern = "^IPM_[A-Z][0-9]",colnames(ALL_2c),value = T)))

#### split fg cols by CROPS and PHI
fg = melt(fg,id.vars = c("LAT"),na.rm = T)
t = colsplit(fg$variable,"_", c("IPM","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## crop 
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## phi

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

fg = fg[,-which(names(fg) %in% "variable")]
fg$crop = t$crop
fg$phi = t$phi
rm(t)
rm(c2)
rm(c1)

# select the tropics
fg$LAT = abs(fg$LAT)
IPM_TROP_2c = dplyr::filter(fg,LAT < 23)

# --- 4c --
fg = dplyr::select(ALL_4c, one_of("LAT",grep(pattern = "^IPM_[A-Z][0-9]",colnames(ALL_4c),value = T)))

#split fg cols by CROPS and PHI
fg = melt(fg,id.vars = c("LAT"),na.rm = T)
t = colsplit(fg$variable,"_", c("IPM","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## crop 
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## phi

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

fg = fg[,-which(names(fg) %in% "variable")]
fg$crop = t$crop
fg$phi = t$phi
rm(t)
rm(c2)
rm(c1)

# select the tropics
fg$LAT = abs(fg$LAT)
IPM_TROP_4c = dplyr::filter(fg,LAT < 23)


# REMOVE OUTLIERS 
IPM_TROP_4c = ddply(IPM_TROP_4c,c("LAT","phi","crop"), function(d){
  limits = median(d$value) + 3*c(-1, 1)*mad(d$value)
  d[(d$value - limits[1])*(limits[2] - d$value) > 0,]
})

IPM_TROP_2c = ddply(IPM_TROP_2c,c("LAT","phi","crop"), function(d){
  limits = median(d$value) + 3*c(-1, 1)*mad(d$value)
  d[(d$value - limits[1])*(limits[2] - d$value) > 0,]
})


# -SUMMARIZE TROPICS (2c and 4c) --
SUM_IPM_TROP_2c = ddply(IPM_TROP_2c,c("phi"),summarise,
              mean= mean(value,na.rm = T),
              median= median(value,na.rm = T),
              min= min(value,na.rm = T),
              max= max(value,na.rm = T))
SUM_IPM_TROP_2c$anomaly = "2C"


SUM_IPM_TROP_4c = ddply(IPM_TROP_4c,c("phi"),summarise,
                         mean= mean(value,na.rm = T),
                         median= median(value,na.rm = T),
                         min= min(value,na.rm = T),
                         max= max(value,na.rm = T))
SUM_IPM_TROP_4c$anomaly = "4C"

# -SUMMARY IPM FOR  THE TROPICS
SUM_IPM_TROPICS = rbind(SUM_IPM_TROP_2c,SUM_IPM_TROP_4c)

rm(fg)
rm(IPM_TROP_4c)
rm(IPM_TROP_2c)
rm(SUM_IPM_TROP_2c)
rm(SUM_IPM_TROP_4c)


# ----5) GET IPM  WITH X Y VALUES ----
# prepare DF 
DAT_2c = ALL_2c[,c("LON","LAT","NAME","region","MET_M2","MET_R2","MET_W2","MET_M3","MET_R3","MET_W3","MET_M4","MET_R4","MET_W4","POP_M2","POP_R2","POP_W2","POP_M3","POP_R3","POP_W3","POP_M4","POP_R4","POP_W4","IPM_M2","IPM_R2","IPM_W2","IPM_M3","IPM_R3","IPM_W3","IPM_M4","IPM_R4","IPM_W4")]

# Change from wide to long
library(reshape2)
DAT_2c = melt(DAT_2c, id.vars=c("LON","LAT","NAME","region"),na.rm=T)

# split cols by fact and crop
t = colsplit(DAT_2c$variable,"_", c("fact","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
t = colsplit(DAT_2c$variable,"_", c("fact","crop"))
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

DAT_2c = DAT_2c[,-which(names(DAT_2c) %in% "variable")]
DAT_2c$crop = t$crop
DAT_2c$phi = t$phi
DAT_2c$fact = t$fact

rm(t)
rm(c2)
rm(c1)

IPM_2c_XY = dplyr::filter(DAT_2c,fact=="IPM")

# save IPM_2c_XY  objects
save(list = "IPM_2c_XY",
     file = file.path(wdrdata,"IPM_2c_XY.RData"))



# --- 7) Map IMP 2c outliers ----
load(file.path(wdrdata,"IPM_2c_XY.RData"))
# summarize IPM outliers
stats.IPM = boxplot.stats(IPM_2c_XY$value)

IPM_2c_OUT = IPM_2c_XY[IPM_2c_XY$value  %in% stats.IPM$out,]

# generate basemap
source(file.path(wdfun,"basemap.r"))

## plot
m1 = ggplot()  + geom_polygon(data=map,aes(x=long, y=lat,group=group), fill="gray", color="gray",size=0.2)
m1 = m1 + geom_raster(data=IPM_2c_OUT,aes(fill=value,x=LON, y=LAT),interpolate = T)
m1 = m1 + scale_fill_gradientn(colours = rev(cols(12)))
m1 = m1 + theme(title=element_blank())
m1 = m1 + ggtitle("IPM Outliers")
m1

# Save plot
ppi = 300
plotname = file.path(wdpng,paste("IPM Outliers",".png",sep = ""))
png(filename=plotname,width=10*ppi, height=8*ppi, res=ppi )
m1
dev.off()

# Density curve
p = ggplot(IPM_2c_OUT, aes(x=value)) + geom_histogram(binwidth = 50)
p = p + ggtitle("IPM Histogram")

ppi = 150
plotname = file.path(wdpng,paste("IPM Histogram",".png",sep = ""))
png(filename=plotname,width=10*ppi, height=8*ppi, res=ppi )
p
dev.off()

# write csv 
write.csv(IPM_2c_OUT, file.path(wdtables,"IPM_OUTLIERS.csv"),row.names=FALSE)
IPM_2c_HIGH0<- subset(IPM_2c_OUT,value>1)
write.csv(IPM_2c_HIGH0, file.path(wdtables,"IPM_HIGH_OUT_2c.csv"),row.names=FALSE)



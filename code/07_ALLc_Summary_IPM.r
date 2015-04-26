#/////////
# IPM SUMMARIES
# evaluate all code and look for the following objects:

# SUM_IPM: IPM Summary across all temp anomalies and phis
# SUM_IPM_REGIONAL: IPM Summary across all temp anomalies, phis and regions
# SUM_IPM_TROPICS: IPM Summary across all temp anomalies and phis within the tropics (absolute LAT < 23)
#////////


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
source(file.path(wdfun,"clean_space.r")) 

library(plyr)
library(dplyr)

# --- EXPLORE OUTLIERS ----
stats.IPM = boxplot.stats(DAT_2c$value[DAT_2c$fact=="IPM"])
length.IPM = length(DAT_2c$value[DAT_2c$fact=="IPM"])
length.IPM.out = length(stats.IPM$out)
percent.out = length.IPM.out / length.IPM
print(percent.out * 100)

boxplot(DAT_2c$value[DAT_2c$fact=="IPM"],main = "IPM Boxplot")

rm(stats.IPM)
rm(length.IPM)
rm(length.IPM.out)
rm(percent.out)


# --- REMOVE IPM OUTLIERS FROM MASTER TABLES ----
DAT_2c = ddply(DAT_2c,c("fact","phi","crop"), function(d){
  limits = median(d$value) + 3*c(-1, 1)*mad(d$value)
  d[(d$value - limits[1])*(limits[2] - d$value) > 0,]
})

DAT_4c = ddply(DAT_4c,c("fact","phi","crop"), function(d){
  limits = median(d$value) + 3*c(-1, 1)*mad(d$value)
  d[(d$value - limits[1])*(limits[2] - d$value) > 0,]
})


#--- SUMMARY BY TEMP ANOM ----
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


#--- REGIONAL IPM SUMMARY  ----
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

#--- TROPICAL IPM SUMMARY  ----

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


# --- REMOVE OUTLIERS ----
IPM_TROP_4c = ddply(IPM_TROP_4c,c("LAT","phi","crop"), function(d){
  limits = median(d$value) + 3*c(-1, 1)*mad(d$value)
  d[(d$value - limits[1])*(limits[2] - d$value) > 0,]
})

IPM_TROP_2c = ddply(IPM_TROP_2c,c("LAT","phi","crop"), function(d){
  limits = median(d$value) + 3*c(-1, 1)*mad(d$value)
  d[(d$value - limits[1])*(limits[2] - d$value) > 0,]
})


# --- SUMMARIZE TROPICS (ALL TEMP) ----
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

# --- SUMMARY IPM FOR  THE TROPICS
SUM_IPM_TROPICS = rbind(SUM_IPM_TROP_2c,SUM_IPM_TROP_4c)

rm(fg)
rm(IPM_TROP_4c)
rm(IPM_TROP_2c)
rm(SUM_IPM_TROP_2c)
rm(SUM_IPM_TROP_4c)

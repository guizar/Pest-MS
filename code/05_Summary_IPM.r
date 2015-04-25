############---
# IPM SUMMARIES
# evaluate all code and look for the following objects:
# SUM_IPM: IPM Summary across all temp anomalies and phis
# SUM_REGION_IPM: IPM Summary across all temp anomalies, phis and regions
# SUM_IPM_TROPICS: IPM Summary across all temp anomalies and phis within the tropics (absolute LAT < 23)
##### ---


# Load data and initial variables ------

rm(list = ls())

wd = "~/R/Pest-MS/"
wdpng = "~/R/Pest-MS/png"
wdtables = "~/R/Pest-MS/tables"
wddata = "~/R/Pest-MS/data/"
wdrdata = "~/R/Pest-MS/RData/"
wdfun = "~/R/Pest-MS/fun"

# Load data
setwd(wdrdata)
load("ALL_2c.RData")
load("ALL_4c.RData")

library(plyr)
library(dplyr)



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
SUM_REGION_IPM = dplyr::filter(SUM,fact=="IPM")





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


#### --- SUMMARIZE TROPICS (ALL TEMP) --- ###
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

# SUMMARY IPM
SUM_IPM_TROPICS = rbind(SUM_IPM_TROP_2c,SUM_IPM_TROP_4c)



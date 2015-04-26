#///////
# Summaries for Supplemental tables
# //////

# ---- SUMMARIZE _COUNTRY_2c 
rm(list=ls())

#---- Load data and initial variables 
wdpng = "~/R/Pest-MS/png"
wdtables = "~/R/Pest-MS/tables"
wdfun = "~/R/Pest-MS/fun"
wdrdata = "~/R/Pest-MS/RData/"

# Load data
setwd(wdrdata)
load("ALL_2c.RData") 

# ---- SUM_DAT_COUNTRY_2c (median, mean, range) ----
library(plyr) 
# summarize by "NAME","crop","phi", "fact" = serves for a general purpose and ggplot
SUM_DAT_COUNTRY_2c <- ddply(DAT_2c, c("fact","NAME","phi","crop"), summarise,
                    median= median(value,na.rm = T),
                    mean= mean(value,na.rm = T),
                    min= min(value,na.rm = T),
                    max= max(value,na.rm = T))

# Save data
write.csv(SUM_DAT_COUNTRY_2c, file.path(wdtables,"SUM_DAT_COUNTRY_2c.csv"),row.names=F)
# save.image(file.path(wdrdata,"ALL_COUNTRY_2c.RData"))



# ---- Table 1 IPM_COUNTRY_2c: Split data by crops + IPM ----
library(dplyr)

# Extract
SUM_IPM_COUNTRY_2c = SUM_DAT_COUNTRY_2c[SUM_DAT_COUNTRY_2c$fact=="IPM",]
# Remove unused
SUM_IPM_COUNTRY_2c = SUM_IPM_COUNTRY_2c[,-which(colnames(SUM_IPM_COUNTRY_2c) %in% c("fact"),arr.ind = T)]
rownames(SUM_IPM_COUNTRY_2c) = NULL

# Covert from long to wide
library(reshape2)
SUM_IPM_COUNTRY_2c = reshape(SUM_IPM_COUNTRY_2c, direction="wide", idvar=c("NAME","crop"), timevar="phi")

#extract crop
IPM_COUNTRY_2c_MAIZE = SUM_IPM_COUNTRY_2c[SUM_IPM_COUNTRY_2c$crop=="M",]
IPM_COUNTRY_2c_RICE = SUM_IPM_COUNTRY_2c[SUM_IPM_COUNTRY_2c$crop=="R",]
IPM_COUNTRY_2c_WHEAT = SUM_IPM_COUNTRY_2c[SUM_IPM_COUNTRY_2c$crop=="W",]

# remove crop col
IPM_COUNTRY_2c_MAIZE = dplyr::select(IPM_COUNTRY_2c_MAIZE,-crop)
IPM_COUNTRY_2c_RICE = dplyr::select(IPM_COUNTRY_2c_RICE,-crop)
IPM_COUNTRY_2c_WHEAT = dplyr::select(IPM_COUNTRY_2c_WHEAT,-crop)

# horizontal re-ordering (columns)
cols = c("NAME","median.2","median.3","median.4","mean.2","mean.3","mean.4","min.2","max.2","min.3","max.3","min.4","max.4")
IPM_COUNTRY_2c_MAIZE = IPM_COUNTRY_2c_MAIZE[,cols]
IPM_COUNTRY_2c_RICE = IPM_COUNTRY_2c_RICE[,cols]
IPM_COUNTRY_2c_WHEAT = IPM_COUNTRY_2c_WHEAT[,cols]

# concatenate min,max columns
IPM_COUNTRY_2c_MAIZE$range.2 =paste(IPM_COUNTRY_2c_MAIZE$min.2,IPM_COUNTRY_2c_MAIZE$max.2,sep=",")
IPM_COUNTRY_2c_MAIZE$range.3 =paste(IPM_COUNTRY_2c_MAIZE$min.3,IPM_COUNTRY_2c_MAIZE$max.3,sep=",")
IPM_COUNTRY_2c_MAIZE$range.4 =paste(IPM_COUNTRY_2c_MAIZE$min.4,IPM_COUNTRY_2c_MAIZE$max.4,sep=",")

IPM_COUNTRY_2c_RICE$range.2 =paste(IPM_COUNTRY_2c_RICE$min.2,IPM_COUNTRY_2c_RICE$max.2,sep=",")
IPM_COUNTRY_2c_RICE$range.3 =paste(IPM_COUNTRY_2c_RICE$min.3,IPM_COUNTRY_2c_RICE$max.3,sep=",")
IPM_COUNTRY_2c_RICE$range.4 =paste(IPM_COUNTRY_2c_RICE$min.4,IPM_COUNTRY_2c_RICE$max.4,sep=",")

IPM_COUNTRY_2c_WHEAT$range.2 =paste(IPM_COUNTRY_2c_WHEAT$min.2,IPM_COUNTRY_2c_WHEAT$max.2,sep=",")
IPM_COUNTRY_2c_WHEAT$range.3 =paste(IPM_COUNTRY_2c_WHEAT$min.3,IPM_COUNTRY_2c_WHEAT$max.3,sep=",")
IPM_COUNTRY_2c_WHEAT$range.4 =paste(IPM_COUNTRY_2c_WHEAT$min.4,IPM_COUNTRY_2c_WHEAT$max.4,sep=",")

# remove min,max columns
rem =  c("min.2","max.2","min.3","max.3","min.4","max.4")
IPM_COUNTRY_2c_MAIZE = dplyr::select(IPM_COUNTRY_2c_MAIZE,-one_of(rem))
IPM_COUNTRY_2c_RICE = dplyr::select(IPM_COUNTRY_2c_RICE,-one_of(rem))
IPM_COUNTRY_2c_WHEAT = dplyr::select(IPM_COUNTRY_2c_WHEAT,-one_of(rem))
rm(rem)


# vertical order by MEDIAN phi 2
IPM_COUNTRY_2c_MAIZE = IPM_COUNTRY_2c_MAIZE[order(IPM_COUNTRY_2c_MAIZE$median.2,decreasing = T),]
IPM_COUNTRY_2c_RICE = IPM_COUNTRY_2c_RICE[order(IPM_COUNTRY_2c_RICE$median.2,decreasing = T),]
IPM_COUNTRY_2c_WHEAT = IPM_COUNTRY_2c_WHEAT[order(IPM_COUNTRY_2c_WHEAT$median.2,decreasing = T),]

write.csv(IPM_COUNTRY_2c_MAIZE,file.path(wdtables,"ST1_IPM_COUNTRY_2c_MAIZE.csv"),row.names=F)
write.csv(IPM_COUNTRY_2c_RICE, file.path(wdtables,"ST1_IPM_COUNTRY_2c_RICE.csv"),row.names=F)
write.csv(IPM_COUNTRY_2c_WHEAT, file.path(wdtables,"ST1_IPM_COUNTRY_2c_WHEAT.csv"),row.names=F)

rm(IPM_COUNTRY_2c_MAIZE)
rm(IPM_COUNTRY_2c_RICE)
rm(IPM_COUNTRY_2c_WHEAT)



# ---- Table 2 TONNES_COUNTRY_2c ----
# Sup 2:  Total yield and projected annual yield loss in 2050, across a full range of insect life-histories 

library(dplyr)

# TOTAL SUM _COUNTRY_2c by "NAME","crop","phi"
SUM_TONNES_COUNTRY_2c <- ddply(TONNES_2c, c("NAME","phi","crop"), summarise,
                    sum= sum(value,na.rm = T))
                    
 
# summarize PRESENT by "NAME","crop",
SUM_TONNES_COUNTRY_PRES <- ddply(TONNES_PRES, c("NAME","crop"), summarise,
                       sum= sum(value,na.rm = T))

# since Guyana is not listed in the 2c table, I will remove it from the PRESENT table as well
SUM_TONNES_COUNTRY_PRES = dplyr::filter(SUM_TONNES_COUNTRY_PRES, NAME != "Guyana")

# Covert from long to wide
library(reshape2)
SUM_TONNES_COUNTRY_2c = reshape(SUM_TONNES_COUNTRY_2c, direction="wide", idvar=c("NAME","crop"), timevar="phi")

#extract crop
TONNES_COUNTRY_2c_MAIZE = SUM_TONNES_COUNTRY_2c[SUM_TONNES_COUNTRY_2c$crop=="M",]
TONNES_COUNTRY_2c_RICE = SUM_TONNES_COUNTRY_2c[SUM_TONNES_COUNTRY_2c$crop=="R",]
TONNES_COUNTRY_2c_WHEAT = SUM_TONNES_COUNTRY_2c[SUM_TONNES_COUNTRY_2c$crop=="W",]

TONNES_COUNTRY_PRES_MAIZE = SUM_TONNES_COUNTRY_PRES[SUM_TONNES_COUNTRY_PRES$crop=="M",]
TONNES_COUNTRY_PRES_RICE = SUM_TONNES_COUNTRY_PRES[SUM_TONNES_COUNTRY_PRES$crop=="R",]
TONNES_COUNTRY_PRES_WHEAT = SUM_TONNES_COUNTRY_PRES[SUM_TONNES_COUNTRY_PRES$crop=="W",]

# remove crop col
TONNES_COUNTRY_2c_MAIZE = TONNES_COUNTRY_2c_MAIZE[,-which(colnames(TONNES_COUNTRY_2c_MAIZE) %in% "crop", arr.ind=T)]
TONNES_COUNTRY_2c_RICE = TONNES_COUNTRY_2c_RICE[,-which(colnames(TONNES_COUNTRY_2c_RICE) %in% "crop", arr.ind=T)]
TONNES_COUNTRY_2c_WHEAT = TONNES_COUNTRY_2c_WHEAT[,-which(colnames(TONNES_COUNTRY_2c_WHEAT) %in% "crop", arr.ind=T)]

# add PRES values _COUNTRY_2c
TONNES_COUNTRY_2c_MAIZE$present =TONNES_COUNTRY_PRES_MAIZE$sum
TONNES_COUNTRY_2c_RICE$present =TONNES_COUNTRY_PRES_RICE$sum
TONNES_COUNTRY_2c_WHEAT$present =TONNES_COUNTRY_PRES_WHEAT$sum

# vertical order by PRES 
TONNES_COUNTRY_2c_MAIZE = TONNES_COUNTRY_2c_MAIZE[order(TONNES_COUNTRY_2c_MAIZE$present,decreasing = T),]
TONNES_COUNTRY_2c_RICE = TONNES_COUNTRY_2c_RICE[order(TONNES_COUNTRY_2c_RICE$present,decreasing = T),]
TONNES_COUNTRY_2c_WHEAT = TONNES_COUNTRY_2c_WHEAT[order(TONNES_COUNTRY_2c_WHEAT$present,decreasing = T),]

# horizontal re-ordering (columns)
cols = c("NAME","present","sum.2","sum.3","sum.4")
TONNES_COUNTRY_2c_MAIZE = TONNES_COUNTRY_2c_MAIZE[,cols]
TONNES_COUNTRY_2c_RICE = TONNES_COUNTRY_2c_RICE[,cols]
TONNES_COUNTRY_2c_WHEAT = TONNES_COUNTRY_2c_WHEAT[,cols]
rm(cols)

write.csv(TONNES_COUNTRY_2c_MAIZE,file.path(wdtables,"ST2_TONNES_COUNTRY_2c_MAIZE.csv"),row.names=F)
write.csv(TONNES_COUNTRY_2c_RICE, file.path(wdtables,"ST2_TONNES_COUNTRY_2c_RICE.csv"),row.names=F)
write.csv(TONNES_COUNTRY_2c_WHEAT, file.path(wdtables,"ST2_TONNES_COUNTRY_2c_WHEAT.csv"),row.names=F)

rm(TONNES_COUNTRY_2c_MAIZE)
rm(TONNES_COUNTRY_2c_RICE)
rm(TONNES_COUNTRY_2c_WHEAT)
rm(TONNES_COUNTRY_PRES_MAIZE)
rm(TONNES_COUNTRY_PRES_RICE)
rm(TONNES_COUNTRY_PRES_WHEAT)



# ---- Table 3 YLPH_COUNTRY_2c ----
# Sup 3: Yield and yield loss (tones per ha) caused by climate change impacts on crop pests,  across a full range of insect life-histories 
library(dplyr)
SUM_YLPH_COUNTRY_2c <- ddply(YLPH_2c, c("NAME","phi","crop"), summarise,
                     median= median(value,na.rm = T),
                     mean= mean(value,na.rm = T),
                     min= min(value,na.rm = T),
                     max= max(value,na.rm = T))


# summarize PRESENT by "NAME","crop",
SUM_YL_COUNTRY_PRES <- ddply(YL_PRES, c("NAME","crop"), summarise,
                     median= median(value,na.rm = T),
                     mean= mean(value,na.rm = T),
                     min= min(value,na.rm = T),
                     max= max(value,na.rm = T))

# Covert from long to wide
library(reshape2)
SUM_YLPH_COUNTRY_2c = reshape(SUM_YLPH_COUNTRY_2c, direction="wide", idvar=c("NAME","crop"), timevar="phi")

#extract crop
YLPH_COUNTRY_2c_MAIZE = SUM_YLPH_COUNTRY_2c[SUM_YLPH_COUNTRY_2c$crop=="M",]
YLPH_COUNTRY_2c_RICE = SUM_YLPH_COUNTRY_2c[SUM_YLPH_COUNTRY_2c$crop=="R",]
YLPH_COUNTRY_2c_WHEAT = SUM_YLPH_COUNTRY_2c[SUM_YLPH_COUNTRY_2c$crop=="W",]

YL_COUNTRY_PRES_MAIZE = SUM_YL_COUNTRY_PRES[SUM_YL_COUNTRY_PRES$crop=="M",]
YL_COUNTRY_PRES_RICE = SUM_YL_COUNTRY_PRES[SUM_YL_COUNTRY_PRES$crop=="R",]
YL_COUNTRY_PRES_WHEAT = SUM_YL_COUNTRY_PRES[SUM_YL_COUNTRY_PRES$crop=="W",]

# remove crop col
YLPH_COUNTRY_2c_MAIZE = dplyr::select(YLPH_COUNTRY_2c_MAIZE,-crop)
YLPH_COUNTRY_2c_RICE = dplyr::select(YLPH_COUNTRY_2c_RICE,-crop)
YLPH_COUNTRY_2c_WHEAT = dplyr::select(YLPH_COUNTRY_2c_WHEAT,-crop)

YL_COUNTRY_PRES_MAIZE = dplyr::select(YL_COUNTRY_PRES_MAIZE,-crop)
YL_COUNTRY_PRES_RICE = dplyr::select(YL_COUNTRY_PRES_RICE,-crop)
YL_COUNTRY_PRES_WHEAT = dplyr::select(YL_COUNTRY_PRES_WHEAT,-crop)

# horizontal re-ordering (columns)
cols = c("NAME","median.2","median.3","median.4","mean.2","mean.3","mean.4","min.2","max.2","min.3","max.3","min.4","max.4")
YLPH_COUNTRY_2c_MAIZE = YLPH_COUNTRY_2c_MAIZE[,cols]
YLPH_COUNTRY_2c_RICE = YLPH_COUNTRY_2c_RICE[,cols]
YLPH_COUNTRY_2c_WHEAT = YLPH_COUNTRY_2c_WHEAT[,cols]
rm(cols)

# concatenate min,max columns
YLPH_COUNTRY_2c_MAIZE$range.2 =paste(YLPH_COUNTRY_2c_MAIZE$min.2,YLPH_COUNTRY_2c_MAIZE$max.2,sep=",")
YLPH_COUNTRY_2c_MAIZE$range.3 =paste(YLPH_COUNTRY_2c_MAIZE$min.3,YLPH_COUNTRY_2c_MAIZE$max.3,sep=",")
YLPH_COUNTRY_2c_MAIZE$range.4 =paste(YLPH_COUNTRY_2c_MAIZE$min.4,YLPH_COUNTRY_2c_MAIZE$max.4,sep=",")

YLPH_COUNTRY_2c_RICE$range.2 =paste(YLPH_COUNTRY_2c_RICE$min.2,YLPH_COUNTRY_2c_RICE$max.2,sep=",")
YLPH_COUNTRY_2c_RICE$range.3 =paste(YLPH_COUNTRY_2c_RICE$min.3,YLPH_COUNTRY_2c_RICE$max.3,sep=",")
YLPH_COUNTRY_2c_RICE$range.4 =paste(YLPH_COUNTRY_2c_RICE$min.4,YLPH_COUNTRY_2c_RICE$max.4,sep=",")

YLPH_COUNTRY_2c_WHEAT$range.2 =paste(YLPH_COUNTRY_2c_WHEAT$min.2,YLPH_COUNTRY_2c_WHEAT$max.2,sep=",")
YLPH_COUNTRY_2c_WHEAT$range.3 =paste(YLPH_COUNTRY_2c_WHEAT$min.3,YLPH_COUNTRY_2c_WHEAT$max.3,sep=",")
YLPH_COUNTRY_2c_WHEAT$range.4 =paste(YLPH_COUNTRY_2c_WHEAT$min.4,YLPH_COUNTRY_2c_WHEAT$max.4,sep=",")

# concatenate min,max columns (PRES)
YL_COUNTRY_PRES_MAIZE$range =paste(YL_COUNTRY_PRES_MAIZE$min,YL_COUNTRY_PRES_MAIZE$max,sep=",")
YL_COUNTRY_PRES_RICE$range =paste(YL_COUNTRY_PRES_RICE$min,YL_COUNTRY_PRES_RICE$max,sep=",")
YL_COUNTRY_PRES_WHEAT$range =paste(YL_COUNTRY_PRES_WHEAT$min,YL_COUNTRY_PRES_WHEAT$max,sep=",")

# remove min,max columns (PRES)
YL_COUNTRY_PRES_MAIZE = dplyr::select(YL_COUNTRY_PRES_MAIZE,-one_of(c("min","max")))
YL_COUNTRY_PRES_RICE = dplyr::select(YL_COUNTRY_PRES_RICE,-one_of(c("min","max")))
YL_COUNTRY_PRES_WHEAT = dplyr::select(YL_COUNTRY_PRES_WHEAT,-one_of(c("min","max")))

# ADD PRESENT VALS TO YLPH
YLPH_COUNTRY_2c_MAIZE = dplyr::right_join(YL_COUNTRY_PRES_MAIZE,YLPH_COUNTRY_2c_MAIZE,by = "NAME")
YLPH_COUNTRY_2c_RICE = dplyr::right_join(YL_COUNTRY_PRES_RICE,YLPH_COUNTRY_2c_RICE,by = "NAME")
YLPH_COUNTRY_2c_WHEAT = dplyr::right_join(YL_COUNTRY_PRES_WHEAT,YLPH_COUNTRY_2c_WHEAT,by = "NAME")

# remove min,max cols (YLPH)
rem =  c("min.2","max.2","min.3","max.3","min.4","max.4")
YLPH_COUNTRY_2c_MAIZE = dplyr::select(YLPH_COUNTRY_2c_MAIZE,-one_of(rem))
YLPH_COUNTRY_2c_RICE = dplyr::select(YLPH_COUNTRY_2c_RICE,-one_of(rem))
YLPH_COUNTRY_2c_WHEAT = dplyr::select(YLPH_COUNTRY_2c_WHEAT,-one_of(rem))
rm(rem)

# arrange by median (present)
YLPH_COUNTRY_2c_MAIZE = dplyr::arrange(YLPH_COUNTRY_2c_MAIZE, desc(median))
YLPH_COUNTRY_2c_RICE = dplyr::arrange(YLPH_COUNTRY_2c_RICE, desc(median))
YLPH_COUNTRY_2c_WHEAT = dplyr::arrange(YLPH_COUNTRY_2c_WHEAT, desc(median))

# write csv
write.csv(YLPH_COUNTRY_2c_MAIZE,file.path(wdtables,"ST3_YLPH_COUNTRY_2c_MAIZE.csv"),row.names=F)
write.csv(YLPH_COUNTRY_2c_RICE, file.path(wdtables,"ST3_YLPH_COUNTRY_2c_RICE.csv"),row.names=F)
write.csv(YLPH_COUNTRY_2c_WHEAT, file.path(wdtables,"ST3_YLPH_COUNTRY_2c_WHEAT.csv"),row.names=F)

# remove unused
rm(YLPH_COUNTRY_2c_MAIZE)
rm(YLPH_COUNTRY_2c_RICE)
rm(YLPH_COUNTRY_2c_WHEAT)
rm(YL_COUNTRY_PRES_MAIZE)
rm(YL_COUNTRY_PRES_RICE)
rm(YL_COUNTRY_PRES_WHEAT)



# ---- Table 4 IYCC_COUNTRY_2c ----
# Table S4: Percent of crop lost due to climate change, across a full range of insect life-histories 
library(dplyr)
library(plyr)

SUM_IYCC_COUNTRY_2c <- ddply(IYCC_2c, c("NAME","phi","crop"), summarise,
                     median= median(value,na.rm = T),
                     mean= mean(value,na.rm = T),
                     min= min(value,na.rm = T),
                     max= max(value,na.rm = T))

# Covert from long to wide
library(reshape2)
SUM_IYCC_COUNTRY_2c = reshape(SUM_IYCC_COUNTRY_2c, direction="wide", idvar=c("NAME","crop"), timevar="phi")

#extract crop
IYCC_COUNTRY_2c_MAIZE = SUM_IYCC_COUNTRY_2c[SUM_IYCC_COUNTRY_2c$crop=="M",]
IYCC_COUNTRY_2c_RICE = SUM_IYCC_COUNTRY_2c[SUM_IYCC_COUNTRY_2c$crop=="R",]
IYCC_COUNTRY_2c_WHEAT = SUM_IYCC_COUNTRY_2c[SUM_IYCC_COUNTRY_2c$crop=="W",]

# remove crop col
IYCC_COUNTRY_2c_MAIZE = dplyr::select(IYCC_COUNTRY_2c_MAIZE,-crop)
IYCC_COUNTRY_2c_RICE = dplyr::select(IYCC_COUNTRY_2c_RICE,-crop)
IYCC_COUNTRY_2c_WHEAT = dplyr::select(IYCC_COUNTRY_2c_WHEAT,-crop)

# horizontal re-ordering (columns)
cols = c("NAME","median.2","median.3","median.4","mean.2","mean.3","mean.4","min.2","max.2","min.3","max.3","min.4","max.4")
IYCC_COUNTRY_2c_MAIZE = IYCC_COUNTRY_2c_MAIZE[,cols]
IYCC_COUNTRY_2c_RICE = IYCC_COUNTRY_2c_RICE[,cols]
IYCC_COUNTRY_2c_WHEAT = IYCC_COUNTRY_2c_WHEAT[,cols]
rm(cols)

# concatenate min,max columns
IYCC_COUNTRY_2c_MAIZE$range.2 =paste(IYCC_COUNTRY_2c_MAIZE$min.2,IYCC_COUNTRY_2c_MAIZE$max.2,sep=",")
IYCC_COUNTRY_2c_MAIZE$range.3 =paste(IYCC_COUNTRY_2c_MAIZE$min.3,IYCC_COUNTRY_2c_MAIZE$max.3,sep=",")
IYCC_COUNTRY_2c_MAIZE$range.4 =paste(IYCC_COUNTRY_2c_MAIZE$min.4,IYCC_COUNTRY_2c_MAIZE$max.4,sep=",")

IYCC_COUNTRY_2c_RICE$range.2 =paste(IYCC_COUNTRY_2c_RICE$min.2,IYCC_COUNTRY_2c_RICE$max.2,sep=",")
IYCC_COUNTRY_2c_RICE$range.3 =paste(IYCC_COUNTRY_2c_RICE$min.3,IYCC_COUNTRY_2c_RICE$max.3,sep=",")
IYCC_COUNTRY_2c_RICE$range.4 =paste(IYCC_COUNTRY_2c_RICE$min.4,IYCC_COUNTRY_2c_RICE$max.4,sep=",")

IYCC_COUNTRY_2c_WHEAT$range.2 =paste(IYCC_COUNTRY_2c_WHEAT$min.2,IYCC_COUNTRY_2c_WHEAT$max.2,sep=",")
IYCC_COUNTRY_2c_WHEAT$range.3 =paste(IYCC_COUNTRY_2c_WHEAT$min.3,IYCC_COUNTRY_2c_WHEAT$max.3,sep=",")
IYCC_COUNTRY_2c_WHEAT$range.4 =paste(IYCC_COUNTRY_2c_WHEAT$min.4,IYCC_COUNTRY_2c_WHEAT$max.4,sep=",")

# remove min,max columns
rem =  c("min.2","max.2","min.3","max.3","min.4","max.4")
IYCC_COUNTRY_2c_MAIZE = dplyr::select(IYCC_COUNTRY_2c_MAIZE,-one_of(rem))
IYCC_COUNTRY_2c_RICE = dplyr::select(IYCC_COUNTRY_2c_RICE,-one_of(rem))
IYCC_COUNTRY_2c_WHEAT = dplyr::select(IYCC_COUNTRY_2c_WHEAT,-one_of(rem))
rm(rem)

# vertical order by MEDIAN phi 2
IYCC_COUNTRY_2c_MAIZE = IYCC_COUNTRY_2c_MAIZE[order(IYCC_COUNTRY_2c_MAIZE$median.2,decreasing = T),]
IYCC_COUNTRY_2c_RICE = IYCC_COUNTRY_2c_RICE[order(IYCC_COUNTRY_2c_RICE$median.2,decreasing = T),]
IYCC_COUNTRY_2c_WHEAT = IYCC_COUNTRY_2c_WHEAT[order(IYCC_COUNTRY_2c_WHEAT$median.2,decreasing = T),]

# write csv
write.csv(IYCC_COUNTRY_2c_MAIZE,file.path(wdtables,"ST4_IYCC_COUNTRY_2c_MAIZE.csv"),row.names=F)
write.csv(IYCC_COUNTRY_2c_RICE, file.path(wdtables,"ST4_IYCC_COUNTRY_2c_RICE.csv"),row.names=F)
write.csv(IYCC_COUNTRY_2c_WHEAT, file.path(wdtables,"ST4_IYCC_COUNTRY_2c_WHEAT.csv"),row.names=F)

# remove unused
rm(IYCC_COUNTRY_2c_MAIZE)
rm(IYCC_COUNTRY_2c_RICE)
rm(IYCC_COUNTRY_2c_WHEAT)

#---- SAVE DATA 2c ----
save.image(file.path(wdrdata,"ALL_COUNTRY_2c.RData"))



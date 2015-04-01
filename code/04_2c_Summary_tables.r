#///////
# Summaries for Suplemental tables
# //////

# ---- SUMMARIZE _2c 
rm(list=ls())

#---- Load data and initial variables 
wdpng = "~/R/Pest-MS/png"
wdtables = "~/R/Pest-MS/tables"
wdfun = "~/R/Pest-MS/functions"
wdrdata = "~/R/Pest-MS/RData/"

# Load data
setwd(wdrdata)
load("ALL_2c.RData") 

# ---- SUM_DAT_2c (median, mean, range) ----
library(plyr) 
# summarize by "region","crop","phi", "fact" = serves for a general purpose and ggplot
SUM_DAT_2c <- ddply(DAT_2c, c("fact","region","phi","crop"), summarise,
                    median= median(value,na.rm = T),
                    mean= mean(value,na.rm = T),
                    min= min(value,na.rm = T),
                    max= max(value,na.rm = T))

# Save data
write.csv(SUM_DAT_2c, file.path(wdtables,"SUM_DAT_2c.csv"),row.names=F)
# save.image(file.path(wdrdata,"ALL_2c.RData"))



# ---- Table 1 IPM_2c: Split data by crops + IPM ----
library(dplyr)

# Extract
SUM_IPM_2c = SUM_DAT_2c[SUM_DAT_2c$fact=="IPM",]
# Remove unused
SUM_IPM_2c = SUM_IPM_2c[,-which(colnames(SUM_IPM_2c) %in% c("fact"),arr.ind = T)]
rownames(SUM_IPM_2c) = NULL

# Covert from long to wide
library(reshape2)
SUM_IPM_2c = reshape(SUM_IPM_2c, direction="wide", idvar=c("region","crop"), timevar="phi")

#extract crop
IPM_2c_MAIZE = SUM_IPM_2c[SUM_IPM_2c$crop=="M",]
IPM_2c_RICE = SUM_IPM_2c[SUM_IPM_2c$crop=="R",]
IPM_2c_WEAT = SUM_IPM_2c[SUM_IPM_2c$crop=="W",]

# remove crop col
IPM_2c_MAIZE = dplyr::select(IPM_2c_MAIZE,-crop)
IPM_2c_RICE = dplyr::select(IPM_2c_RICE,-crop)
IPM_2c_WEAT = dplyr::select(IPM_2c_WEAT,-crop)

# horizontal re-ordering (columns)
cols = c("region","median.2","median.3","median.4","mean.2","mean.3","mean.4","min.2","max.2","min.3","max.3","min.4","max.4")
IPM_2c_MAIZE = IPM_2c_MAIZE[,cols]
IPM_2c_RICE = IPM_2c_RICE[,cols]
IPM_2c_WEAT = IPM_2c_WEAT[,cols]

# concatenate min,max columns
IPM_2c_MAIZE$range.2 =paste(IPM_2c_MAIZE$min.2,IPM_2c_MAIZE$max.2,sep=",")
IPM_2c_MAIZE$range.3 =paste(IPM_2c_MAIZE$min.3,IPM_2c_MAIZE$max.3,sep=",")
IPM_2c_MAIZE$range.4 =paste(IPM_2c_MAIZE$min.4,IPM_2c_MAIZE$max.4,sep=",")

IPM_2c_RICE$range.2 =paste(IPM_2c_RICE$min.2,IPM_2c_RICE$max.2,sep=",")
IPM_2c_RICE$range.3 =paste(IPM_2c_RICE$min.3,IPM_2c_RICE$max.3,sep=",")
IPM_2c_RICE$range.4 =paste(IPM_2c_RICE$min.4,IPM_2c_RICE$max.4,sep=",")

IPM_2c_WEAT$range.2 =paste(IPM_2c_WEAT$min.2,IPM_2c_WEAT$max.2,sep=",")
IPM_2c_WEAT$range.3 =paste(IPM_2c_WEAT$min.3,IPM_2c_WEAT$max.3,sep=",")
IPM_2c_WEAT$range.4 =paste(IPM_2c_WEAT$min.4,IPM_2c_WEAT$max.4,sep=",")

# remove min,max columns
rem =  c("min.2","max.2","min.3","max.3","min.4","max.4")
IPM_2c_MAIZE = dplyr::select(IPM_2c_MAIZE,-one_of(rem))
IPM_2c_RICE = dplyr::select(IPM_2c_RICE,-one_of(rem))
IPM_2c_WEAT = dplyr::select(IPM_2c_WEAT,-one_of(rem))
rm(rem)


# vertical order by MEDIAN phi 2
IPM_2c_MAIZE = IPM_2c_MAIZE[order(IPM_2c_MAIZE$median.2,decreasing = T),]
IPM_2c_RICE = IPM_2c_RICE[order(IPM_2c_RICE$median.2,decreasing = T),]
IPM_2c_WEAT = IPM_2c_WEAT[order(IPM_2c_WEAT$median.2,decreasing = T),]

write.csv(IPM_2c_MAIZE,file.path(wdtables,"ST1_IPM_2c_MAIZE.csv"),row.names=F)
write.csv(IPM_2c_RICE, file.path(wdtables,"ST1_IPM_2c_RICE.csv"),row.names=F)
write.csv(IPM_2c_WEAT, file.path(wdtables,"ST1_IPM_2c_WEAT.csv"),row.names=F)

rm(IPM_2c_MAIZE)
rm(IPM_2c_RICE)
rm(IPM_2c_WEAT)



# ---- Table 2 TONNES_2c ----
# Sup 2:  Total yield and projected annual yield loss in 2050, across a full range of insect life-histories 

library(dplyr)

# TOTAL SUM _2c by "region","crop","phi"
SUM_TONNES_2c <- ddply(TONNES_2c, c("region","phi","crop"), summarise,
                    sum= sum(value,na.rm = T))
                    

# summarize PRESENT by "region","crop",
SUM_TONNES_PRES <- ddply(TONNES_PRES, c("region","crop"), summarise,
                       sum= sum(value,na.rm = T))

# Covert from long to wide
library(reshape2)
SUM_TONNES_2c = reshape(SUM_TONNES_2c, direction="wide", idvar=c("region","crop"), timevar="phi")

#extract crop
TONNES_2c_MAIZE = SUM_TONNES_2c[SUM_TONNES_2c$crop=="M",]
TONNES_2c_RICE = SUM_TONNES_2c[SUM_TONNES_2c$crop=="R",]
TONNES_2c_WEAT = SUM_TONNES_2c[SUM_TONNES_2c$crop=="W",]

TONNES_PRES_MAIZE = SUM_TONNES_PRES[SUM_TONNES_PRES$crop=="M",]
TONNES_PRES_RICE = SUM_TONNES_PRES[SUM_TONNES_PRES$crop=="R",]
TONNES_PRES_WEAT = SUM_TONNES_PRES[SUM_TONNES_PRES$crop=="W",]

# remove crop col
TONNES_2c_MAIZE = TONNES_2c_MAIZE[,-which(colnames(TONNES_2c_MAIZE) %in% "crop", arr.ind=T)]
TONNES_2c_RICE = TONNES_2c_RICE[,-which(colnames(TONNES_2c_RICE) %in% "crop", arr.ind=T)]
TONNES_2c_WEAT = TONNES_2c_WEAT[,-which(colnames(TONNES_2c_WEAT) %in% "crop", arr.ind=T)]

# add PRES values _2c
TONNES_2c_MAIZE$present =TONNES_PRES_MAIZE$sum
TONNES_2c_RICE$present =TONNES_PRES_RICE$sum
TONNES_2c_WEAT$present =TONNES_PRES_WEAT$sum

# vertical order by PRES 
TONNES_2c_MAIZE = TONNES_2c_MAIZE[order(TONNES_2c_MAIZE$present,decreasing = T),]
TONNES_2c_RICE = TONNES_2c_RICE[order(TONNES_2c_RICE$present,decreasing = T),]
TONNES_2c_WEAT = TONNES_2c_WEAT[order(TONNES_2c_WEAT$present,decreasing = T),]

# horizontal re-ordering (columns)
cols = c("region","present","sum.2","sum.3","sum.4")
TONNES_2c_MAIZE = TONNES_2c_MAIZE[,cols]
TONNES_2c_MAIZE = TONNES_2c_MAIZE[,cols]
TONNES_2c_MAIZE = TONNES_2c_MAIZE[,cols]
rm(cols)

write.csv(TONNES_2c_MAIZE,file.path(wdtables,"ST2_TONNES_2c_MAIZE.csv"),row.names=F)
write.csv(TONNES_2c_RICE, file.path(wdtables,"ST2_TONNES_2c_RICE.csv"),row.names=F)
write.csv(TONNES_2c_WEAT, file.path(wdtables,"ST2_TONNES_2c_WEAT.csv"),row.names=F)

rm(TONNES_2c_MAIZE)
rm(TONNES_2c_RICE)
rm(TONNES_2c_WEAT)
rm(TONNES_PRES_MAIZE)
rm(TONNES_PRES_RICE)
rm(TONNES_PRES_WEAT)



# ---- Table 3 YLPH_2c ----
# Sup 3: Yield and yield loss (tones per ha) caused by climate change impacts on crop pests,  across a full range of insect life-histories 
library(dplyr)
SUM_YLPH_2c <- ddply(YLPH_2c, c("region","phi","crop"), summarise,
                     median= median(value,na.rm = T),
                     mean= mean(value,na.rm = T),
                     min= min(value,na.rm = T),
                     max= max(value,na.rm = T))


# summarize PRESENT by "region","crop",
SUM_YL_PRES <- ddply(YL_PRES, c("region","crop"), summarise,
                     median= median(value,na.rm = T),
                     mean= mean(value,na.rm = T),
                     min= min(value,na.rm = T),
                     max= max(value,na.rm = T))

# Covert from long to wide
library(reshape2)
SUM_YLPH_2c = reshape(SUM_YLPH_2c, direction="wide", idvar=c("region","crop"), timevar="phi")

#extract crop
YLPH_2c_MAIZE = SUM_YLPH_2c[SUM_YLPH_2c$crop=="M",]
YLPH_2c_RICE = SUM_YLPH_2c[SUM_YLPH_2c$crop=="R",]
YLPH_2c_WEAT = SUM_YLPH_2c[SUM_YLPH_2c$crop=="W",]

YL_PRES_MAIZE = SUM_YL_PRES[SUM_YL_PRES$crop=="M",]
YL_PRES_RICE = SUM_YL_PRES[SUM_YL_PRES$crop=="R",]
YL_PRES_WEAT = SUM_YL_PRES[SUM_YL_PRES$crop=="W",]

# remove crop col
YLPH_2c_MAIZE = dplyr::select(YLPH_2c_MAIZE,-crop)
YLPH_2c_RICE = dplyr::select(YLPH_2c_RICE,-crop)
YLPH_2c_WEAT = dplyr::select(YLPH_2c_WEAT,-crop)

YL_PRES_MAIZE = dplyr::select(YL_PRES_MAIZE,-crop)
YL_PRES_RICE = dplyr::select(YL_PRES_RICE,-crop)
YL_PRES_WEAT = dplyr::select(YL_PRES_WEAT,-crop)

# horizontal re-ordering (columns)
cols = c("region","median.2","median.3","median.4","mean.2","mean.3","mean.4","min.2","max.2","min.3","max.3","min.4","max.4")
YLPH_2c_MAIZE = YLPH_2c_MAIZE[,cols]
YLPH_2c_RICE = YLPH_2c_RICE[,cols]
YLPH_2c_WEAT = YLPH_2c_WEAT[,cols]
rm(cols)

# concatenate min,max columns
YLPH_2c_MAIZE$range.2 =paste(YLPH_2c_MAIZE$min.2,YLPH_2c_MAIZE$max.2,sep=",")
YLPH_2c_MAIZE$range.3 =paste(YLPH_2c_MAIZE$min.3,YLPH_2c_MAIZE$max.3,sep=",")
YLPH_2c_MAIZE$range.4 =paste(YLPH_2c_MAIZE$min.4,YLPH_2c_MAIZE$max.4,sep=",")

YLPH_2c_RICE$range.2 =paste(YLPH_2c_RICE$min.2,YLPH_2c_RICE$max.2,sep=",")
YLPH_2c_RICE$range.3 =paste(YLPH_2c_RICE$min.3,YLPH_2c_RICE$max.3,sep=",")
YLPH_2c_RICE$range.4 =paste(YLPH_2c_RICE$min.4,YLPH_2c_RICE$max.4,sep=",")

YLPH_2c_WEAT$range.2 =paste(YLPH_2c_WEAT$min.2,YLPH_2c_WEAT$max.2,sep=",")
YLPH_2c_WEAT$range.3 =paste(YLPH_2c_WEAT$min.3,YLPH_2c_WEAT$max.3,sep=",")
YLPH_2c_WEAT$range.4 =paste(YLPH_2c_WEAT$min.4,YLPH_2c_WEAT$max.4,sep=",")

# concatenate min,max columns (PRES)
YL_PRES_MAIZE$range =paste(YL_PRES_MAIZE$min,YL_PRES_MAIZE$max,sep=",")
YL_PRES_RICE$range =paste(YL_PRES_RICE$min,YL_PRES_RICE$max,sep=",")
YL_PRES_WEAT$range =paste(YL_PRES_WEAT$min,YL_PRES_WEAT$max,sep=",")

# remove min,max columns (PRES)
YL_PRES_MAIZE = dplyr::select(YL_PRES_MAIZE,-one_of(c("min","max")))
YL_PRES_RICE = dplyr::select(YL_PRES_RICE,-one_of(c("min","max")))
YL_PRES_WEAT = dplyr::select(YL_PRES_WEAT,-one_of(c("min","max")))

# ADD PRESENT VALS TO YLPH
YLPH_2c_MAIZE = dplyr::right_join(YL_PRES_MAIZE,YLPH_2c_MAIZE,by = "region")
YLPH_2c_RICE = dplyr::right_join(YL_PRES_RICE,YLPH_2c_RICE,by = "region")
YLPH_2c_WEAT = dplyr::right_join(YL_PRES_WEAT,YLPH_2c_WEAT,by = "region")

# remove min,max cols (YLPH)
rem =  c("min.2","max.2","min.3","max.3","min.4","max.4")
YLPH_2c_MAIZE = dplyr::select(YLPH_2c_MAIZE,-one_of(rem))
YLPH_2c_RICE = dplyr::select(YLPH_2c_RICE,-one_of(rem))
YLPH_2c_WEAT = dplyr::select(YLPH_2c_WEAT,-one_of(rem))
rm(rem)

# arrange by median (present)
YLPH_2c_MAIZE = dplyr::arrange(YLPH_2c_MAIZE, desc(median))
YLPH_2c_RICE = dplyr::arrange(YLPH_2c_RICE, desc(median))
YLPH_2c_WEAT = dplyr::arrange(YLPH_2c_WEAT, desc(median))

# write csv
write.csv(YLPH_2c_MAIZE,file.path(wdtables,"ST3_YLPH_2c_MAIZE.csv"),row.names=F)
write.csv(YLPH_2c_RICE, file.path(wdtables,"ST3_YLPH_2c_RICE.csv"),row.names=F)
write.csv(YLPH_2c_WEAT, file.path(wdtables,"ST3_YLPH_2c_WEAT.csv"),row.names=F)

# remove unused
rm(YLPH_2c_MAIZE)
rm(YLPH_2c_RICE)
rm(YLPH_2c_WEAT)
rm(YL_PRES_MAIZE)
rm(YL_PRES_RICE)
rm(YL_PRES_WEAT)



# ---- Table 4 IYCC_2c ----
# Table S4: Percent of crop lost due to climate change, across a full range of insect life-histories 
library(dplyr)
library(plyr)

SUM_IYCC_2c <- ddply(IYCC_2c, c("region","phi","crop"), summarise,
                     median= median(value,na.rm = T),
                     mean= mean(value,na.rm = T),
                     min= min(value,na.rm = T),
                     max= max(value,na.rm = T))

# Covert from long to wide
library(reshape2)
SUM_IYCC_2c = reshape(SUM_IYCC_2c, direction="wide", idvar=c("region","crop"), timevar="phi")

#extract crop
IYCC_2c_MAIZE = SUM_IYCC_2c[SUM_IYCC_2c$crop=="M",]
IYCC_2c_RICE = SUM_IYCC_2c[SUM_IYCC_2c$crop=="R",]
IYCC_2c_WEAT = SUM_IYCC_2c[SUM_IYCC_2c$crop=="W",]

# remove crop col
IYCC_2c_MAIZE = dplyr::select(IYCC_2c_MAIZE,-crop)
IYCC_2c_RICE = dplyr::select(IYCC_2c_RICE,-crop)
IYCC_2c_WEAT = dplyr::select(IYCC_2c_WEAT,-crop)

# horizontal re-ordering (columns)
cols = c("region","median.2","median.3","median.4","mean.2","mean.3","mean.4","min.2","max.2","min.3","max.3","min.4","max.4")
IYCC_2c_MAIZE = IYCC_2c_MAIZE[,cols]
IYCC_2c_RICE = IYCC_2c_RICE[,cols]
IYCC_2c_WEAT = IYCC_2c_WEAT[,cols]
rm(cols)

# concatenate min,max columns
IYCC_2c_MAIZE$range.2 =paste(IYCC_2c_MAIZE$min.2,IYCC_2c_MAIZE$max.2,sep=",")
IYCC_2c_MAIZE$range.3 =paste(IYCC_2c_MAIZE$min.3,IYCC_2c_MAIZE$max.3,sep=",")
IYCC_2c_MAIZE$range.4 =paste(IYCC_2c_MAIZE$min.4,IYCC_2c_MAIZE$max.4,sep=",")

IYCC_2c_RICE$range.2 =paste(IYCC_2c_RICE$min.2,IYCC_2c_RICE$max.2,sep=",")
IYCC_2c_RICE$range.3 =paste(IYCC_2c_RICE$min.3,IYCC_2c_RICE$max.3,sep=",")
IYCC_2c_RICE$range.4 =paste(IYCC_2c_RICE$min.4,IYCC_2c_RICE$max.4,sep=",")

IYCC_2c_WEAT$range.2 =paste(IYCC_2c_WEAT$min.2,IYCC_2c_WEAT$max.2,sep=",")
IYCC_2c_WEAT$range.3 =paste(IYCC_2c_WEAT$min.3,IYCC_2c_WEAT$max.3,sep=",")
IYCC_2c_WEAT$range.4 =paste(IYCC_2c_WEAT$min.4,IYCC_2c_WEAT$max.4,sep=",")

# remove min,max columns
rem =  c("min.2","max.2","min.3","max.3","min.4","max.4")
IYCC_2c_MAIZE = dplyr::select(IYCC_2c_MAIZE,-one_of(rem))
IYCC_2c_RICE = dplyr::select(IYCC_2c_RICE,-one_of(rem))
IYCC_2c_WEAT = dplyr::select(IYCC_2c_WEAT,-one_of(rem))
rm(rem)

# vertical order by MEDIAN phi 2
IYCC_2c_MAIZE = IYCC_2c_MAIZE[order(IYCC_2c_MAIZE$median.2,decreasing = T),]
IYCC_2c_RICE = IYCC_2c_RICE[order(IYCC_2c_RICE$median.2,decreasing = T),]
IYCC_2c_WEAT = IYCC_2c_WEAT[order(IYCC_2c_WEAT$median.2,decreasing = T),]

# write csv
write.csv(IYCC_2c_MAIZE,file.path(wdtables,"ST4_IYCC_2c_MAIZE.csv"),row.names=F)
write.csv(IYCC_2c_RICE, file.path(wdtables,"ST4_IYCC_2c_RICE.csv"),row.names=F)
write.csv(IYCC_2c_WEAT, file.path(wdtables,"ST4_IYCC_2c_WEAT.csv"),row.names=F)

# remove unused
rm(IYCC_2c_MAIZE)
rm(IYCC_2c_RICE)
rm(IYCC_2c_WEAT)

#---- SAVE DATA 2c ----
save.image(file.path(wdrdata,"ALL_2c.RData"))



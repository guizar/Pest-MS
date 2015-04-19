#///////
# Summaries for Supplemental tables
# //////

# ---- SUMMARIZE _REGION_4c 
rm(list=ls())

#---- Load data and initial variables 
wdpng = "~/R/Pest-MS/png"
wdtables = "~/R/Pest-MS/tables"
wdfun = "~/R/Pest-MS/fun"
wdrdata = "~/R/Pest-MS/RData/"

# Load data
setwd(wdrdata)
load("ALL_4c.RData") 

# ---- SUM_DAT_REGION_4c (median, mean, range) ----
library(plyr) 
# summarize by "region","crop","phi", "fact" = serves for a general purpose and ggplot
SUM_DAT_REGION_4c <- ddply(DAT_4c, c("fact","region","phi","crop"), summarise,
                    median= median(value,na.rm = T),
                    mean= mean(value,na.rm = T),
                    min= min(value,na.rm = T),
                    max= max(value,na.rm = T))

# Save data
write.csv(SUM_DAT_REGION_4c, file.path(wdtables,"SUM_DAT_REGION_4c.csv"),row.names=F)
# save.image(file.path(wdrdata,"ALL_REGION_4c.RData"))



# ---- Table 1 IPM_REGION_4c: Split data by crops + IPM ----
library(dplyr)

# Extract
SUM_IPM_REGION_4c = SUM_DAT_REGION_4c[SUM_DAT_REGION_4c$fact=="IPM",]
# Remove unused
SUM_IPM_REGION_4c = SUM_IPM_REGION_4c[,-which(colnames(SUM_IPM_REGION_4c) %in% c("fact"),arr.ind = T)]
rownames(SUM_IPM_REGION_4c) = NULL

# Covert from long to wide
library(reshape2)
SUM_IPM_REGION_4c = reshape(SUM_IPM_REGION_4c, direction="wide", idvar=c("region","crop"), timevar="phi")

#extract crop
IPM_REGION_4c_MAIZE = SUM_IPM_REGION_4c[SUM_IPM_REGION_4c$crop=="M",]
IPM_REGION_4c_RICE = SUM_IPM_REGION_4c[SUM_IPM_REGION_4c$crop=="R",]
IPM_REGION_4c_WEAT = SUM_IPM_REGION_4c[SUM_IPM_REGION_4c$crop=="W",]

# remove crop col
IPM_REGION_4c_MAIZE = dplyr::select(IPM_REGION_4c_MAIZE,-crop)
IPM_REGION_4c_RICE = dplyr::select(IPM_REGION_4c_RICE,-crop)
IPM_REGION_4c_WEAT = dplyr::select(IPM_REGION_4c_WEAT,-crop)

# horizontal re-ordering (columns)
cols = c("region","median.2","median.3","median.4","mean.2","mean.3","mean.4","min.2","max.2","min.3","max.3","min.4","max.4")
IPM_REGION_4c_MAIZE = IPM_REGION_4c_MAIZE[,cols]
IPM_REGION_4c_RICE = IPM_REGION_4c_RICE[,cols]
IPM_REGION_4c_WEAT = IPM_REGION_4c_WEAT[,cols]

# concatenate min,max columns
IPM_REGION_4c_MAIZE$range.2 =paste(IPM_REGION_4c_MAIZE$min.2,IPM_REGION_4c_MAIZE$max.2,sep=",")
IPM_REGION_4c_MAIZE$range.3 =paste(IPM_REGION_4c_MAIZE$min.3,IPM_REGION_4c_MAIZE$max.3,sep=",")
IPM_REGION_4c_MAIZE$range.4 =paste(IPM_REGION_4c_MAIZE$min.4,IPM_REGION_4c_MAIZE$max.4,sep=",")

IPM_REGION_4c_RICE$range.2 =paste(IPM_REGION_4c_RICE$min.2,IPM_REGION_4c_RICE$max.2,sep=",")
IPM_REGION_4c_RICE$range.3 =paste(IPM_REGION_4c_RICE$min.3,IPM_REGION_4c_RICE$max.3,sep=",")
IPM_REGION_4c_RICE$range.4 =paste(IPM_REGION_4c_RICE$min.4,IPM_REGION_4c_RICE$max.4,sep=",")

IPM_REGION_4c_WEAT$range.2 =paste(IPM_REGION_4c_WEAT$min.2,IPM_REGION_4c_WEAT$max.2,sep=",")
IPM_REGION_4c_WEAT$range.3 =paste(IPM_REGION_4c_WEAT$min.3,IPM_REGION_4c_WEAT$max.3,sep=",")
IPM_REGION_4c_WEAT$range.4 =paste(IPM_REGION_4c_WEAT$min.4,IPM_REGION_4c_WEAT$max.4,sep=",")

# remove min,max columns
rem =  c("min.2","max.2","min.3","max.3","min.4","max.4")
IPM_REGION_4c_MAIZE = dplyr::select(IPM_REGION_4c_MAIZE,-one_of(rem))
IPM_REGION_4c_RICE = dplyr::select(IPM_REGION_4c_RICE,-one_of(rem))
IPM_REGION_4c_WEAT = dplyr::select(IPM_REGION_4c_WEAT,-one_of(rem))
rm(rem)


# vertical order by MEDIAN phi 2
IPM_REGION_4c_MAIZE = IPM_REGION_4c_MAIZE[order(IPM_REGION_4c_MAIZE$median.2,decreasing = T),]
IPM_REGION_4c_RICE = IPM_REGION_4c_RICE[order(IPM_REGION_4c_RICE$median.2,decreasing = T),]
IPM_REGION_4c_WEAT = IPM_REGION_4c_WEAT[order(IPM_REGION_4c_WEAT$median.2,decreasing = T),]

write.csv(IPM_REGION_4c_MAIZE,file.path(wdtables,"ST1_IPM_REGION_4c_MAIZE.csv"),row.names=F)
write.csv(IPM_REGION_4c_RICE, file.path(wdtables,"ST1_IPM_REGION_4c_RICE.csv"),row.names=F)
write.csv(IPM_REGION_4c_WEAT, file.path(wdtables,"ST1_IPM_REGION_4c_WEAT.csv"),row.names=F)

rm(IPM_REGION_4c_MAIZE)
rm(IPM_REGION_4c_RICE)
rm(IPM_REGION_4c_WEAT)



# ---- Table 2 TONNES_REGION_4c ----
# Sup 2:  Total yield and projected annual yield loss in 2050, across a full range of insect life-histories 

library(dplyr)

# TOTAL SUM _REGION_4c by "region","crop","phi"
SUM_TONNES_REGION_4c <- ddply(TONNES_4c, c("region","phi","crop"), summarise,
                    sum= sum(value,na.rm = T))
                    

# summarize PRESENT by "region","crop",
SUM_TONNES_REGION_PRES <- ddply(TONNES_PRES, c("region","crop"), summarise,
                       sum= sum(value,na.rm = T))

# Covert from long to wide
library(reshape2)
SUM_TONNES_REGION_4c = reshape(SUM_TONNES_REGION_4c, direction="wide", idvar=c("region","crop"), timevar="phi")

#extract crop
TONNES_REGION_4c_MAIZE = SUM_TONNES_REGION_4c[SUM_TONNES_REGION_4c$crop=="M",]
TONNES_REGION_4c_RICE = SUM_TONNES_REGION_4c[SUM_TONNES_REGION_4c$crop=="R",]
TONNES_REGION_4c_WEAT = SUM_TONNES_REGION_4c[SUM_TONNES_REGION_4c$crop=="W",]

TONNES_REGION_PRES_MAIZE = SUM_TONNES_REGION_PRES[SUM_TONNES_REGION_PRES$crop=="M",]
TONNES_REGION_PRES_RICE = SUM_TONNES_REGION_PRES[SUM_TONNES_REGION_PRES$crop=="R",]
TONNES_REGION_PRES_WEAT = SUM_TONNES_REGION_PRES[SUM_TONNES_REGION_PRES$crop=="W",]

# remove crop col
TONNES_REGION_4c_MAIZE = TONNES_REGION_4c_MAIZE[,-which(colnames(TONNES_REGION_4c_MAIZE) %in% "crop", arr.ind=T)]
TONNES_REGION_4c_RICE = TONNES_REGION_4c_RICE[,-which(colnames(TONNES_REGION_4c_RICE) %in% "crop", arr.ind=T)]
TONNES_REGION_4c_WEAT = TONNES_REGION_4c_WEAT[,-which(colnames(TONNES_REGION_4c_WEAT) %in% "crop", arr.ind=T)]

# add PRES values _REGION_4c
TONNES_REGION_4c_MAIZE$present =TONNES_REGION_PRES_MAIZE$sum
TONNES_REGION_4c_RICE$present =TONNES_REGION_PRES_RICE$sum
TONNES_REGION_4c_WEAT$present =TONNES_REGION_PRES_WEAT$sum

# vertical order by PRES 
TONNES_REGION_4c_MAIZE = TONNES_REGION_4c_MAIZE[order(TONNES_REGION_4c_MAIZE$present,decreasing = T),]
TONNES_REGION_4c_RICE = TONNES_REGION_4c_RICE[order(TONNES_REGION_4c_RICE$present,decreasing = T),]
TONNES_REGION_4c_WEAT = TONNES_REGION_4c_WEAT[order(TONNES_REGION_4c_WEAT$present,decreasing = T),]

# horizontal re-ordering (columns)
cols = c("region","present","sum.2","sum.3","sum.4")
TONNES_REGION_4c_MAIZE = TONNES_REGION_4c_MAIZE[,cols]
TONNES_REGION_4c_RICE = TONNES_REGION_4c_RICE[,cols]
TONNES_REGION_4c_WEAT = TONNES_REGION_4c_WEAT[,cols]
rm(cols)

write.csv(TONNES_REGION_4c_MAIZE,file.path(wdtables,"ST2_TONNES_REGION_4c_MAIZE.csv"),row.names=F)
write.csv(TONNES_REGION_4c_RICE, file.path(wdtables,"ST2_TONNES_REGION_4c_RICE.csv"),row.names=F)
write.csv(TONNES_REGION_4c_WEAT, file.path(wdtables,"ST2_TONNES_REGION_4c_WEAT.csv"),row.names=F)

rm(TONNES_REGION_4c_MAIZE)
rm(TONNES_REGION_4c_RICE)
rm(TONNES_REGION_4c_WEAT)
rm(TONNES_REGION_PRES_MAIZE)
rm(TONNES_REGION_PRES_RICE)
rm(TONNES_REGION_PRES_WEAT)



# ---- Table 3 YLPH_REGION_4c ----
# Sup 3: Yield and yield loss (tones per ha) caused by climate change impacts on crop pests,  across a full range of insect life-histories 
library(dplyr)
SUM_YLPH_REGION_4c <- ddply(YLPH_4c, c("region","phi","crop"), summarise,
                     median= median(value,na.rm = T),
                     mean= mean(value,na.rm = T),
                     min= min(value,na.rm = T),
                     max= max(value,na.rm = T))


# summarize PRESENT by "region","crop",
SUM_YL_REGION_PRES <- ddply(YL_PRES, c("region","crop"), summarise,
                     median= median(value,na.rm = T),
                     mean= mean(value,na.rm = T),
                     min= min(value,na.rm = T),
                     max= max(value,na.rm = T))

# Covert from long to wide
library(reshape2)
SUM_YLPH_REGION_4c = reshape(SUM_YLPH_REGION_4c, direction="wide", idvar=c("region","crop"), timevar="phi")

#extract crop
YLPH_REGION_4c_MAIZE = SUM_YLPH_REGION_4c[SUM_YLPH_REGION_4c$crop=="M",]
YLPH_REGION_4c_RICE = SUM_YLPH_REGION_4c[SUM_YLPH_REGION_4c$crop=="R",]
YLPH_REGION_4c_WEAT = SUM_YLPH_REGION_4c[SUM_YLPH_REGION_4c$crop=="W",]

YL_REGION_PRES_MAIZE = SUM_YL_REGION_PRES[SUM_YL_REGION_PRES$crop=="M",]
YL_REGION_PRES_RICE = SUM_YL_REGION_PRES[SUM_YL_REGION_PRES$crop=="R",]
YL_REGION_PRES_WEAT = SUM_YL_REGION_PRES[SUM_YL_REGION_PRES$crop=="W",]

# remove crop col
YLPH_REGION_4c_MAIZE = dplyr::select(YLPH_REGION_4c_MAIZE,-crop)
YLPH_REGION_4c_RICE = dplyr::select(YLPH_REGION_4c_RICE,-crop)
YLPH_REGION_4c_WEAT = dplyr::select(YLPH_REGION_4c_WEAT,-crop)

YL_REGION_PRES_MAIZE = dplyr::select(YL_REGION_PRES_MAIZE,-crop)
YL_REGION_PRES_RICE = dplyr::select(YL_REGION_PRES_RICE,-crop)
YL_REGION_PRES_WEAT = dplyr::select(YL_REGION_PRES_WEAT,-crop)

# horizontal re-ordering (columns)
cols = c("region","median.2","median.3","median.4","mean.2","mean.3","mean.4","min.2","max.2","min.3","max.3","min.4","max.4")
YLPH_REGION_4c_MAIZE = YLPH_REGION_4c_MAIZE[,cols]
YLPH_REGION_4c_RICE = YLPH_REGION_4c_RICE[,cols]
YLPH_REGION_4c_WEAT = YLPH_REGION_4c_WEAT[,cols]
rm(cols)

# concatenate min,max columns
YLPH_REGION_4c_MAIZE$range.2 =paste(YLPH_REGION_4c_MAIZE$min.2,YLPH_REGION_4c_MAIZE$max.2,sep=",")
YLPH_REGION_4c_MAIZE$range.3 =paste(YLPH_REGION_4c_MAIZE$min.3,YLPH_REGION_4c_MAIZE$max.3,sep=",")
YLPH_REGION_4c_MAIZE$range.4 =paste(YLPH_REGION_4c_MAIZE$min.4,YLPH_REGION_4c_MAIZE$max.4,sep=",")

YLPH_REGION_4c_RICE$range.2 =paste(YLPH_REGION_4c_RICE$min.2,YLPH_REGION_4c_RICE$max.2,sep=",")
YLPH_REGION_4c_RICE$range.3 =paste(YLPH_REGION_4c_RICE$min.3,YLPH_REGION_4c_RICE$max.3,sep=",")
YLPH_REGION_4c_RICE$range.4 =paste(YLPH_REGION_4c_RICE$min.4,YLPH_REGION_4c_RICE$max.4,sep=",")

YLPH_REGION_4c_WEAT$range.2 =paste(YLPH_REGION_4c_WEAT$min.2,YLPH_REGION_4c_WEAT$max.2,sep=",")
YLPH_REGION_4c_WEAT$range.3 =paste(YLPH_REGION_4c_WEAT$min.3,YLPH_REGION_4c_WEAT$max.3,sep=",")
YLPH_REGION_4c_WEAT$range.4 =paste(YLPH_REGION_4c_WEAT$min.4,YLPH_REGION_4c_WEAT$max.4,sep=",")

# concatenate min,max columns (PRES)
YL_REGION_PRES_MAIZE$range =paste(YL_REGION_PRES_MAIZE$min,YL_REGION_PRES_MAIZE$max,sep=",")
YL_REGION_PRES_RICE$range =paste(YL_REGION_PRES_RICE$min,YL_REGION_PRES_RICE$max,sep=",")
YL_REGION_PRES_WEAT$range =paste(YL_REGION_PRES_WEAT$min,YL_REGION_PRES_WEAT$max,sep=",")

# remove min,max columns (PRES)
YL_REGION_PRES_MAIZE = dplyr::select(YL_REGION_PRES_MAIZE,-one_of(c("min","max")))
YL_REGION_PRES_RICE = dplyr::select(YL_REGION_PRES_RICE,-one_of(c("min","max")))
YL_REGION_PRES_WEAT = dplyr::select(YL_REGION_PRES_WEAT,-one_of(c("min","max")))

# ADD PRESENT VALS TO YLPH
YLPH_REGION_4c_MAIZE = dplyr::right_join(YL_REGION_PRES_MAIZE,YLPH_REGION_4c_MAIZE,by = "region")
YLPH_REGION_4c_RICE = dplyr::right_join(YL_REGION_PRES_RICE,YLPH_REGION_4c_RICE,by = "region")
YLPH_REGION_4c_WEAT = dplyr::right_join(YL_REGION_PRES_WEAT,YLPH_REGION_4c_WEAT,by = "region")

# remove min,max cols (YLPH)
rem =  c("min.2","max.2","min.3","max.3","min.4","max.4")
YLPH_REGION_4c_MAIZE = dplyr::select(YLPH_REGION_4c_MAIZE,-one_of(rem))
YLPH_REGION_4c_RICE = dplyr::select(YLPH_REGION_4c_RICE,-one_of(rem))
YLPH_REGION_4c_WEAT = dplyr::select(YLPH_REGION_4c_WEAT,-one_of(rem))
rm(rem)

# arrange by median (present)
YLPH_REGION_4c_MAIZE = dplyr::arrange(YLPH_REGION_4c_MAIZE, desc(median))
YLPH_REGION_4c_RICE = dplyr::arrange(YLPH_REGION_4c_RICE, desc(median))
YLPH_REGION_4c_WEAT = dplyr::arrange(YLPH_REGION_4c_WEAT, desc(median))

# write csv
write.csv(YLPH_REGION_4c_MAIZE,file.path(wdtables,"ST3_YLPH_REGION_4c_MAIZE.csv"),row.names=F)
write.csv(YLPH_REGION_4c_RICE, file.path(wdtables,"ST3_YLPH_REGION_4c_RICE.csv"),row.names=F)
write.csv(YLPH_REGION_4c_WEAT, file.path(wdtables,"ST3_YLPH_REGION_4c_WEAT.csv"),row.names=F)

# remove unused
rm(YLPH_REGION_4c_MAIZE)
rm(YLPH_REGION_4c_RICE)
rm(YLPH_REGION_4c_WEAT)
rm(YL_REGION_PRES_MAIZE)
rm(YL_REGION_PRES_RICE)
rm(YL_REGION_PRES_WEAT)



# ---- Table 4 IYCC_REGION_4c ----
# Table S4: Percent of crop lost due to climate change, across a full range of insect life-histories 
library(dplyr)
library(plyr)

SUM_IYCC_REGION_4c <- ddply(IYCC_4c, c("region","phi","crop"), summarise,
                     median= median(value,na.rm = T),
                     mean= mean(value,na.rm = T),
                     min= min(value,na.rm = T),
                     max= max(value,na.rm = T))

# Covert from long to wide
library(reshape2)
SUM_IYCC_REGION_4c = reshape(SUM_IYCC_REGION_4c, direction="wide", idvar=c("region","crop"), timevar="phi")

#extract crop
IYCC_REGION_4c_MAIZE = SUM_IYCC_REGION_4c[SUM_IYCC_REGION_4c$crop=="M",]
IYCC_REGION_4c_RICE = SUM_IYCC_REGION_4c[SUM_IYCC_REGION_4c$crop=="R",]
IYCC_REGION_4c_WEAT = SUM_IYCC_REGION_4c[SUM_IYCC_REGION_4c$crop=="W",]

# remove crop col
IYCC_REGION_4c_MAIZE = dplyr::select(IYCC_REGION_4c_MAIZE,-crop)
IYCC_REGION_4c_RICE = dplyr::select(IYCC_REGION_4c_RICE,-crop)
IYCC_REGION_4c_WEAT = dplyr::select(IYCC_REGION_4c_WEAT,-crop)

# horizontal re-ordering (columns)
cols = c("region","median.2","median.3","median.4","mean.2","mean.3","mean.4","min.2","max.2","min.3","max.3","min.4","max.4")
IYCC_REGION_4c_MAIZE = IYCC_REGION_4c_MAIZE[,cols]
IYCC_REGION_4c_RICE = IYCC_REGION_4c_RICE[,cols]
IYCC_REGION_4c_WEAT = IYCC_REGION_4c_WEAT[,cols]
rm(cols)

# concatenate min,max columns
IYCC_REGION_4c_MAIZE$range.2 =paste(IYCC_REGION_4c_MAIZE$min.2,IYCC_REGION_4c_MAIZE$max.2,sep=",")
IYCC_REGION_4c_MAIZE$range.3 =paste(IYCC_REGION_4c_MAIZE$min.3,IYCC_REGION_4c_MAIZE$max.3,sep=",")
IYCC_REGION_4c_MAIZE$range.4 =paste(IYCC_REGION_4c_MAIZE$min.4,IYCC_REGION_4c_MAIZE$max.4,sep=",")

IYCC_REGION_4c_RICE$range.2 =paste(IYCC_REGION_4c_RICE$min.2,IYCC_REGION_4c_RICE$max.2,sep=",")
IYCC_REGION_4c_RICE$range.3 =paste(IYCC_REGION_4c_RICE$min.3,IYCC_REGION_4c_RICE$max.3,sep=",")
IYCC_REGION_4c_RICE$range.4 =paste(IYCC_REGION_4c_RICE$min.4,IYCC_REGION_4c_RICE$max.4,sep=",")

IYCC_REGION_4c_WEAT$range.2 =paste(IYCC_REGION_4c_WEAT$min.2,IYCC_REGION_4c_WEAT$max.2,sep=",")
IYCC_REGION_4c_WEAT$range.3 =paste(IYCC_REGION_4c_WEAT$min.3,IYCC_REGION_4c_WEAT$max.3,sep=",")
IYCC_REGION_4c_WEAT$range.4 =paste(IYCC_REGION_4c_WEAT$min.4,IYCC_REGION_4c_WEAT$max.4,sep=",")

# remove min,max columns
rem =  c("min.2","max.2","min.3","max.3","min.4","max.4")
IYCC_REGION_4c_MAIZE = dplyr::select(IYCC_REGION_4c_MAIZE,-one_of(rem))
IYCC_REGION_4c_RICE = dplyr::select(IYCC_REGION_4c_RICE,-one_of(rem))
IYCC_REGION_4c_WEAT = dplyr::select(IYCC_REGION_4c_WEAT,-one_of(rem))
rm(rem)

# vertical order by MEDIAN phi 2
IYCC_REGION_4c_MAIZE = IYCC_REGION_4c_MAIZE[order(IYCC_REGION_4c_MAIZE$median.2,decreasing = T),]
IYCC_REGION_4c_RICE = IYCC_REGION_4c_RICE[order(IYCC_REGION_4c_RICE$median.2,decreasing = T),]
IYCC_REGION_4c_WEAT = IYCC_REGION_4c_WEAT[order(IYCC_REGION_4c_WEAT$median.2,decreasing = T),]

# write csv
write.csv(IYCC_REGION_4c_MAIZE,file.path(wdtables,"ST4_IYCC_REGION_4c_MAIZE.csv"),row.names=F)
write.csv(IYCC_REGION_4c_RICE, file.path(wdtables,"ST4_IYCC_REGION_4c_RICE.csv"),row.names=F)
write.csv(IYCC_REGION_4c_WEAT, file.path(wdtables,"ST4_IYCC_REGION_4c_WEAT.csv"),row.names=F)

# remove unused
rm(IYCC_REGION_4c_MAIZE)
rm(IYCC_REGION_4c_RICE)
rm(IYCC_REGION_4c_WEAT)

#---- SAVE DATA 4c ----
save.image(file.path(wdrdata,"ALL_REGION_4c.RData"))



#///////
# Prepare data frames for both graphics and tables
# //////

# ---- ALL_2c (DAT_2c) ----
# Used in graphics and in the Sup table: 1
wddata = "~/R/Pest-MS/data/"
wdrdata = "~/R/Pest-MS/RData/"
setwd(wdrdata)

# load Rdata 
load("ALL_2c.RData") 

# prepare DF 
DAT_2c = ALL_2c[,c("region","MET_M2","MET_R2","MET_W2","MET_M3","MET_R3","MET_W3","MET_M4","MET_R4","MET_W4","POP_M2","POP_R2","POP_W2","POP_M3","POP_R3","POP_W3","POP_M4","POP_R4","POP_W4","IPM_M2","IPM_R2","IPM_W2","IPM_M3","IPM_R3","IPM_W3","IPM_M4","IPM_R4","IPM_W4")]

# Change from wide to long
library(reshape2)
DAT_2c = melt(DAT_2c, id.vars="region",na.rm=T)

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

# save data 
# write.csv(DAT_2c, file.path(wddata,"DAT_2c.csv"),row.names=FALSE)
# save.image(file.path(wdrdata,"ALL_2c.RData"))


# ---- TONNES_2c ----
# Sup 2:  Total yield and projected annual yield loss in 2050, across a full range of insect life-histories 

# YIELD (tonnes) TOTAL PER CELL (YLD_TOT_x)
TONNES_PRES = ALL_2c[,c("region","YLD_TOT_M","YLD_TOT_R","YLD_TOT_W")]

# TOTAL TONNES LOST PER CELL DUE TO CLIMATE INDUCED CHANGES IN INSECT PESTS (CROP LOSS PER CELL- CL2050_xy)
TONNES_2c = ALL_2c[,c("region","CL2050_M2","CL2050_M3","CL2050_M4","CL2050_R2","CL2050_R3","CL2050_R4","CL2050_W2","CL2050_W3","CL2050_W4")]

# Change from wide to long
library(reshape2)
TONNES_PRES = melt(TONNES_PRES, id.vars="region",na.rm=T)
TONNES_2c = melt(TONNES_2c, id.vars="region",na.rm=T)

# split PRES cols by CROP
t1 = colsplit(TONNES_PRES$variable,"_", c("YIELD_TOT","crop"))
t1 = colsplit(t1$crop,"_", c("YIELD_TOT","crop"))
TONNES_PRES$crop = t1$crop
TONNES_PRES = TONNES_PRES[,-which(colnames(TONNES_PRES) %in% "variable")]
rm(t1)

# split TONNES_2c cols by CROP
t = colsplit(TONNES_2c$variable,"_", c("fact","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

TONNES_2c = TONNES_2c[,-which(names(TONNES_2c) %in% "variable")]
TONNES_2c$crop = t$crop
TONNES_2c$phi = t$phi

rm(t)
rm(c2)
rm(c1)

# save data 
# write.csv(TONNES_2c, file.path(wddata,"TONNES_2c.csv"),row.names=FALSE)
# write.csv(TONNES_PRES, file.path(wddata,"TONNES_PRES.csv"),row.names=FALSE)


# ---- YLPH_2c----
# Sup 3: Yield and yield loss (tones per ha) caused by climate change impacts on crop pests,  across a full range of insect life-histories 

library(plyr)
# Select relevant columns
YLPH_2c = dplyr::select(ALL_2c,region,IPM_M2,IPM_M3,IPM_M4,IPM_R2,IPM_R3,IPM_R4,IPM_W2,IPM_W3,IPM_W4,CLF_M,CY_M,CA_M,CLF_R,CY_R,CA_R,CLF_W,CY_W,CA_W)

# Compute yield (tonnes) lost per ha
YLPH_2c <- dplyr::mutate(YLPH_2c,
               YLPH_M2=(IPM_M2*CLF_M*(CY_M/CA_M)), 
               YLPH_M3=(IPM_M3*CLF_M*(CY_M/CA_M)),
               YLPH_M4=(IPM_M4*CLF_M*(CY_M/CA_M)),
               YLPH_R2=(IPM_R2*CLF_R*(CY_R/CA_R)), 
               YLPH_R3=(IPM_R3*CLF_R*(CY_R/CA_R)),
               YLPH_R4=(IPM_R4*CLF_R*(CY_R/CA_R)),
               YLPH_W2=(IPM_W2*CLF_W*(CY_W/CA_W)), 
               YLPH_W3=(IPM_W3*CLF_W*(CY_W/CA_W)),
               YLPH_W4=(IPM_W4*CLF_W*(CY_W/CA_W)))

# Remove unused cols
YLPH_2c = dplyr::select(YLPH_2c,-c(IPM_M2,IPM_M3,IPM_M4,IPM_R2,IPM_R3,IPM_R4,IPM_W2,IPM_W3,IPM_W4,CLF_M,CY_M,CA_M,CLF_R,CY_R,CA_R,CLF_W,CY_W,CA_W))

# PRESENT yield
YL_PRES = dplyr::select(ALL_2c,region,CY_M,CY_R,CY_W)


# Produce long version of the tables
# Change from wide to long
library(reshape2)
YL_PRES = melt(YL_PRES, id.vars="region",na.rm=T)
YLPH_2c = melt(YLPH_2c, id.vars="region",na.rm=T)

# split PRES cols by CROP
t1 = colsplit(YL_PRES$variable,"_", c("CY","crop"))
YL_PRES$crop = t1$crop
YL_PRES = YL_PRES[,-which(colnames(YL_PRES) %in% "variable")]
rm(t1)

# split YLPH_2c cols by CROP
t = colsplit(YLPH_2c$variable,"_", c("YLPH","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

YLPH_2c = YLPH_2c[,-which(names(YLPH_2c) %in% "variable")]
YLPH_2c$crop = t$crop
YLPH_2c$phi = t$phi
rm(t)
rm(c2)
rm(c1)

# save data 
# write.csv(YLPH_2c, file.path(wddata,"YLPH_2c.csv"),row.names=FALSE)
# write.csv(YL_PRES, file.path(wddata,"YL_PRES.csv"),row.names=FALSE)


# ---- IYCC_2c----
# Table S4: Percent of crop lost due to climate change, across a full range of insect life-histories

library(plyr)
# Select relevant columns
IYCC_2c = dplyr::select(ALL_2c,region,IYCC_M2,IYCC_M3,IYCC_M4,IYCC_R2,IYCC_R3,IYCC_R4,IYCC_W2,IYCC_W3,IYCC_W4)

# Produce long version of the tables
# Change from wide to long
library(reshape2)
IYCC_2c = melt(IYCC_2c, id.vars="region",na.rm=T)

# split IYCC cols by CROP
t1 = colsplit(IYCC_2c$variable,"_", c("IYCC","crop"))
c2 = colsplit(t1$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
c1 = colsplit(t1$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t1$crop = c1$crop
t1$phi = c1$phi

IYCC_2c = IYCC_2c[,-which(names(IYCC_2c) %in% "variable")]
IYCC_2c$crop = t1$crop
IYCC_2c$phi = t1$phi
rm(t1)
rm(c2)
rm(c1)

# save data 
# write.csv(YLPH_2c, file.path(wddata,"YLPH_2c.csv"),row.names=FALSE)
# write.csv(YL_PRES, file.path(wddata,"YL_PRES.csv"),row.names=FALSE)


#---- SAVE DATA 2c ----
save.image(file.path(wdrdata,"ALL_2c.RData"))






# ---- ALL_4c (DAT_4c) ----
rm(list=ls())

# Used in graphics and in the Sup table: 1
wddata = "~/R/Pest-MS/data/"
wdrdata = "~/R/Pest-MS/RData/"
setwd(wdrdata)

# load Rdata 
load("ALL_4c.RData") 

# prepare DF 
DAT_4c = ALL_4c[,c("region","MET_M2","MET_R2","MET_W2","MET_M3","MET_R3","MET_W3","MET_M4","MET_R4","MET_W4","POP_M2","POP_R2","POP_W2","POP_M3","POP_R3","POP_W3","POP_M4","POP_R4","POP_W4","IPM_M2","IPM_R2","IPM_W2","IPM_M3","IPM_R3","IPM_W3","IPM_M4","IPM_R4","IPM_W4")]

# Change from wide to long
library(reshape2)
DAT_4c = melt(DAT_4c, id.vars="region",na.rm=T)

# split cols by fact and crop
t = colsplit(DAT_4c$variable,"_", c("fact","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
t = colsplit(DAT_4c$variable,"_", c("fact","crop"))
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

DAT_4c = DAT_4c[,-which(names(DAT_4c) %in% "variable")]
DAT_4c$crop = t$crop
DAT_4c$phi = t$phi
DAT_4c$fact = t$fact

rm(t)
rm(c2)
rm(c1)

# save data 
# write.csv(DAT_4c, file.path(wddata,"DAT_4c.csv"),row.names=FALSE)
# save.image(file.path(wdrdata,"ALL_4c.RData"))


# ---- TONNES_4c ----
# Sup 2:  Total yield and projected annual yield loss in 2050, across a full range of insect life-histories 

# YIELD (tonnes) TOTAL PER CELL (YLD_TOT_x)
TONNES_PRES = ALL_4c[,c("region","YLD_TOT_M","YLD_TOT_R","YLD_TOT_W")]

# TOTAL TONNES LOST PER CELL DUE TO CLIMATE INDUCED CHANGES IN INSECT PESTS (CROP LOSS PER CELL- CL2050_xy)
TONNES_4c = ALL_4c[,c("region","CL2050_M2","CL2050_M3","CL2050_M4","CL2050_R2","CL2050_R3","CL2050_R4","CL2050_W2","CL2050_W3","CL2050_W4")]

# Change from wide to long
library(reshape2)
TONNES_PRES = melt(TONNES_PRES, id.vars="region",na.rm=T)
TONNES_4c = melt(TONNES_4c, id.vars="region",na.rm=T)

# split PRES cols by CROP
t1 = colsplit(TONNES_PRES$variable,"_", c("YIELD_TOT","crop"))
t1 = colsplit(t1$crop,"_", c("YIELD_TOT","crop"))
TONNES_PRES$crop = t1$crop
TONNES_PRES = TONNES_PRES[,-which(colnames(TONNES_PRES) %in% "variable")]
rm(t1)

# split TONNES_4c cols by CROP
t = colsplit(TONNES_4c$variable,"_", c("fact","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

TONNES_4c = TONNES_4c[,-which(names(TONNES_4c) %in% "variable")]
TONNES_4c$crop = t$crop
TONNES_4c$phi = t$phi

rm(t)
rm(c2)
rm(c1)

# save data 
# write.csv(TONNES_4c, file.path(wddata,"TONNES_4c.csv"),row.names=FALSE)
# write.csv(TONNES_PRES, file.path(wddata,"TONNES_PRES.csv"),row.names=FALSE)


# ---- YLPH_4c----
# Sup 3: Yield and yield loss (tones per ha) caused by climate change impacts on crop pests,  across a full range of insect life-histories 

library(plyr)
# Select relevant columns
YLPH_4c = dplyr::select(ALL_4c,region,IPM_M2,IPM_M3,IPM_M4,IPM_R2,IPM_R3,IPM_R4,IPM_W2,IPM_W3,IPM_W4,CLF_M,CY_M,CA_M,CLF_R,CY_R,CA_R,CLF_W,CY_W,CA_W)

# Compute yield (tonnes) lost per ha
YLPH_4c <- dplyr::mutate(YLPH_4c,
                         YLPH_M2=(IPM_M2*CLF_M*(CY_M/CA_M)), 
                         YLPH_M3=(IPM_M3*CLF_M*(CY_M/CA_M)),
                         YLPH_M4=(IPM_M4*CLF_M*(CY_M/CA_M)),
                         YLPH_R2=(IPM_R2*CLF_R*(CY_R/CA_R)), 
                         YLPH_R3=(IPM_R3*CLF_R*(CY_R/CA_R)),
                         YLPH_R4=(IPM_R4*CLF_R*(CY_R/CA_R)),
                         YLPH_W2=(IPM_W2*CLF_W*(CY_W/CA_W)), 
                         YLPH_W3=(IPM_W3*CLF_W*(CY_W/CA_W)),
                         YLPH_W4=(IPM_W4*CLF_W*(CY_W/CA_W)))

# Remove unused cols
YLPH_4c = dplyr::select(YLPH_4c,-c(IPM_M2,IPM_M3,IPM_M4,IPM_R2,IPM_R3,IPM_R4,IPM_W2,IPM_W3,IPM_W4,CLF_M,CY_M,CA_M,CLF_R,CY_R,CA_R,CLF_W,CY_W,CA_W))

# PRESENT yield
YL_PRES = dplyr::select(ALL_4c,region,CY_M,CY_R,CY_W)


# Produce long version of the tables
# Change from wide to long
library(reshape2)
YL_PRES = melt(YL_PRES, id.vars="region",na.rm=T)
YLPH_4c = melt(YLPH_4c, id.vars="region",na.rm=T)

# split PRES cols by CROP
t1 = colsplit(YL_PRES$variable,"_", c("CY","crop"))
YL_PRES$crop = t1$crop
YL_PRES = YL_PRES[,-which(colnames(YL_PRES) %in% "variable")]
rm(t1)

# split YLPH_4c cols by CROP
t = colsplit(YLPH_4c$variable,"_", c("YLPH","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

YLPH_4c = YLPH_4c[,-which(names(YLPH_4c) %in% "variable")]
YLPH_4c$crop = t$crop
YLPH_4c$phi = t$phi
rm(t)
rm(c2)
rm(c1)

# save data 
# write.csv(YLPH_4c, file.path(wddata,"YLPH_4c.csv"),row.names=FALSE)
# write.csv(YL_PRES, file.path(wddata,"YL_PRES.csv"),row.names=FALSE)


# ---- IYCC_4c----
# Table S4: Percent of crop lost due to climate change, across a full range of insect life-histories

library(plyr)
# Select relevant columns
IYCC_4c = dplyr::select(ALL_4c,region,IYCC_M2,IYCC_M3,IYCC_M4,IYCC_R2,IYCC_R3,IYCC_R4,IYCC_W2,IYCC_W3,IYCC_W4)

# Produce long version of the tables
# Change from wide to long
library(reshape2)
IYCC_4c = melt(IYCC_4c, id.vars="region",na.rm=T)

# split IYCC cols by CROP
t1 = colsplit(IYCC_4c$variable,"_", c("IYCC","crop"))
c2 = colsplit(t1$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
c1 = colsplit(t1$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t1$crop = c1$crop
t1$phi = c1$phi

IYCC_4c = IYCC_4c[,-which(names(IYCC_4c) %in% "variable")]
IYCC_4c$crop = t1$crop
IYCC_4c$phi = t1$phi
rm(t1)
rm(c2)
rm(c1)

# save data 
# write.csv(YLPH_4c, file.path(wddata,"YLPH_4c.csv"),row.names=FALSE)
# write.csv(YL_PRES, file.path(wddata,"YL_PRES.csv"),row.names=FALSE)


#---- SAVE DATA 4c ----
save.image(file.path(wdrdata,"ALL_4c.RData"))

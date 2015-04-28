library(doBy)

#///////
# Summary stats for Manuscript
# //////

rm(list = ls())

# ---- ALL_2c (DAT_2csum) ----
# Used in graphics and in the Sup table: 1
wdtables = "~/R/Pest-MS/tables/"
wdrdata = "~/R/Pest-MS/RData/"

# load Rdata 
load(file.path(wdrdata,"ALL_2c.RData")) 

# prepare DF 
DAT_2csum = ALL_2c[,c("NAME","region","MET_M2","MET_R2","MET_W2","MET_M3","MET_R3","MET_W3","MET_M4","MET_R4","MET_W4","POP_M2","POP_R2","POP_W2","POP_M3","POP_R3","POP_W3","POP_M4","POP_R4","POP_W4","IPM_M2","IPM_R2","IPM_W2","IPM_M3","IPM_R3","IPM_W3","IPM_M4","IPM_R4","IPM_W4","IPM_AVG2","IPM_AVG3","IPM_AVG4")]

# Change from wide to long
library(reshape2)
DAT_2csum = melt(DAT_2csum, id.vars=c("NAME","region"),na.rm=T)

# split cols by fact and crop
t = colsplit(DAT_2csum$variable,"_", c("fact","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
t = colsplit(DAT_2csum$variable,"_", c("fact","crop"))
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

DAT_2csum = DAT_2csum[,-which(names(DAT_2csum) %in% "variable")]
DAT_2csum$crop = t$crop
DAT_2csum$phi = t$phi
DAT_2csum$fact = t$fact

rm(t)
rm(c2)
rm(c1)

# ---- ALL_4c (DAT_4csum) ----
# Used in graphics and in the Sup table: 1
wdtables = "~/R/Pest-MS/tables/"
wdrdata = "~/R/Pest-MS/RData/"


# load Rdata 
load(file.path(wdrdata,"ALL_4c.RData")) 

# prepare DF 
DAT_4csum = ALL_4c[,c("NAME","region","MET_M2","MET_R2","MET_W2","MET_M3","MET_R3","MET_W3","MET_M4","MET_R4","MET_W4","POP_M2","POP_R2","POP_W2","POP_M3","POP_R3","POP_W3","POP_M4","POP_R4","POP_W4","IPM_M2","IPM_R2","IPM_W2","IPM_M3","IPM_R3","IPM_W3","IPM_M4","IPM_R4","IPM_W4","IPM_AVG2","IPM_AVG3","IPM_AVG4")]

# Change from wide to long
library(reshape2)
DAT_4csum = melt(DAT_4csum, id.vars=c("NAME","region"),na.rm=T)

# split cols by fact and crop
t = colsplit(DAT_4csum$variable,"_", c("fact","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
t = colsplit(DAT_4csum$variable,"_", c("fact","crop"))
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

DAT_4csum = DAT_4csum[,-which(names(DAT_4csum) %in% "variable")]
DAT_4csum$crop = t$crop
DAT_4csum$phi = t$phi
DAT_4csum$fact = t$fact

rm(t)
rm(c2)
rm(c1)

# Load data
load(file.path(wdrdata,"DAT_2csum.RData"))
load(file.path(wdrdata,"DAT_4csum.RData"))




#median values (global) - 2 degree
summaryBy(value ~ crop + phi + fact, data = DAT_2csum, FUN = function(x) { c(med = median(x), mean = mean(x), sd = sd(x)) } )
#median values (global) - 4 degree
summaryBy(value ~ crop + phi + fact, data = DAT_4csum, FUN = function(x) { c(med = median(x), mean = mean(x), sd = sd(x)) } )


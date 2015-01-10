#///////
# Group data by regions and split columns, preserving "phi","crop" and "factor"
# //////

# Load ALL CSV
wddata = "~/R/Pest-MS/data"
setwd(wddata)
ALL<-read.csv("ALL.csv", header=TRUE,na.strings = c("#VALUE!", "#N/A", "N/A", "NA", ""))

summary = ALL[,c("region","MET_M2","MET_R2","MET_W2","MET_M3","MET_R3","MET_W3","MET_M4","MET_R4","MET_W4","POP_M2","POP_R2","POP_W2","POP_M3","POP_R3","POP_W3","POP_M4","POP_R4","POP_W4","IPM_M2","IPM_R2","IPM_W2","IPM_M3","IPM_R3","IPM_W3","IPM_M4","IPM_R4","IPM_W4")]

# Change from wide to long
library(reshape2)
summary = melt(summary, id.vars="region")

# split cols
t = colsplit(summary$variable,"_", c("fact","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi"))
t = colsplit(summary$variable,"_", c("fact","crop"))
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi"))

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

summary = summary[,-which(names(summary) %in% "variable")]
summary$crop = t$crop
summary$phi = t$phi
summary$fact = t$fact

rm(t)
rm(c2)
rm(c1)

# save data 
# write.csv(summary, file.path(wddata,"summary.csv"))
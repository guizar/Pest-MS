library(doBy)
library(plyr)
library(dplyr)
library(tidyr)

#///////
# Summary stats for Manuscript
# //////

rm(list = ls())

# ---- Load Master datasets
wdtables = "~/R/Pest-MS/tables/"
wdrdata = "~/R/Pest-MS/RData/"

# load Rdata 
load(file.path(wdrdata,"ALL_2c.RData")) 
load(file.path(wdrdata,"ALL_4c.RData")) 

# clean space - (leave ALL and DAT objects only)
source(file.path(wdfun,"clean_space.r")) 

# --- 1) [CACHED - go to step 3] Summaries ALL_2c and ALL_4c ----
# DAT_2csum = ALL_2c[,c("NAME","region","MET_M2","MET_R2","MET_W2","MET_M3","MET_R3","MET_W3","MET_M4","MET_R4","MET_W4","POP_M2","POP_R2","POP_W2","POP_M3","POP_R3","POP_W3","POP_M4","POP_R4","POP_W4","IPM_M2","IPM_R2","IPM_W2","IPM_M3","IPM_R3","IPM_W3","IPM_M4","IPM_R4","IPM_W4","IPM_AVG2","IPM_AVG3","IPM_AVG4")]
# 
# # Change from wide to long
# library(reshape2)
# DAT_2csum = melt(DAT_2csum, id.vars=c("NAME","region"),na.rm=T)
# 
# # split cols by fact and crop
# t = colsplit(DAT_2csum$variable,"_", c("fact","crop"))
# c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
# t = colsplit(DAT_2csum$variable,"_", c("fact","crop"))
# c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values
# 
# # join new columns and delete unused 
# c1$crop = c2$crop
# t$crop = c1$crop
# t$phi = c1$phi
# 
# DAT_2csum = DAT_2csum[,-which(names(DAT_2csum) %in% "variable")]
# DAT_2csum$crop = t$crop
# DAT_2csum$phi = t$phi
# DAT_2csum$fact = t$fact
# 
# rm(t)
# rm(c2)
# rm(c1)
# 
# DAT_4csum = ALL_4c[,c("NAME","region","MET_M2","MET_R2","MET_W2","MET_M3","MET_R3","MET_W3","MET_M4","MET_R4","MET_W4","POP_M2","POP_R2","POP_W2","POP_M3","POP_R3","POP_W3","POP_M4","POP_R4","POP_W4","IPM_M2","IPM_R2","IPM_W2","IPM_M3","IPM_R3","IPM_W3","IPM_M4","IPM_R4","IPM_W4","IPM_AVG2","IPM_AVG3","IPM_AVG4")]
# 
# # Change from wide to long
# library(reshape2)
# DAT_4csum = melt(DAT_4csum, id.vars=c("NAME","region"),na.rm=T)
# 
# # split cols by fact and crop
# t = colsplit(DAT_4csum$variable,"_", c("fact","crop"))
# c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
# t = colsplit(DAT_4csum$variable,"_", c("fact","crop"))
# c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values
# 
# # join new columns and delete unused 
# c1$crop = c2$crop
# t$crop = c1$crop
# t$phi = c1$phi
# 
# DAT_4csum = DAT_4csum[,-which(names(DAT_4csum) %in% "variable")]
# DAT_4csum$crop = t$crop
# DAT_4csum$phi = t$phi
# DAT_4csum$fact = t$fact
# 
# rm(t)
# rm(c2)
# rm(c1)

# --- 2) [CACHED - go to step 3] median values (global) 2 and 4 degrees  ----
# 2 degree
# DAT_2csum = summaryBy(value ~ crop + phi + fact, data = DAT_2csum, FUN = function(x) { c(med = median(x), mean = mean(x), sd = sd(x)) } )
 
# # 4 degree
# DAT_4csum = summaryBy(value ~ crop + phi + fact, data = DAT_4csum, FUN = function(x) { c(med = median(x), mean = mean(x), sd = sd(x)) } )

# # save DAT sum objects
# save(list = ls()[which(ls() %in% c(ls(pattern ="^DAT_4csum"),
#                                    ls(pattern ="^DAT_2csum")))],
#      file = file.path(wdrdata,"DAT_summaries.RData"))

# --- 3) Load DAT sum objects----
load(file.path(wdrdata,"DAT_summaries.RData")) 

# --- 4) National yield quantiles per crop
load(file.path(wdrdata,"ALL_2c.RData")) 

# clean space
rm(list = ls()[
  -which(ls() %in% c(ls(pattern ="^ALL_"),
                     ls(pattern ="^DAT_"),
                     ls(pattern ="^TONNES_PRES")
                     ))])

TONNES = ddply(TONNES_PRES, c("NAME","crop"), summarise,
                  tonnes= sum(value,na.rm = T))

# Total yield by country
YLD_COUNTRY_STATS =  TONNES %>% 
  group_by(NAME) %>%
  transmute(total = sum(tonnes)) %>%  
  group_by(NAME,total) %>%
  filter(row_number(NAME) == 1) %>%
  arrange(desc(total))
  
# Country percentiles 
t = YLD_COUNTRY_STATS[,2] %>% 
  mutate_each(funs(percent_rank(.)))

YLD_COUNTRY_STATS$percentile = t$total

rm(t)

# Percentiles per crop
YLD_PERCENTILES_CROP =  TONNES %>% 
  spread(crop,tonnes)

t = YLD_PERCENTILES_CROP[,2:4] %>% 
  mutate_each(funs(percent_rank(.)))

YLD_PERCENTILES_CROP$M = t$M
YLD_PERCENTILES_CROP$R = t$R
YLD_PERCENTILES_CROP$W = t$W

# arrange by crops
YLD_PERCENTILES_M = YLD_PERCENTILES_CROP[,c("NAME","M")]
YLD_PERCENTILES_M = YLD_PERCENTILES_M[order(YLD_PERCENTILES_M$M, decreasing = T),]

YLD_PERCENTILES_R = YLD_PERCENTILES_CROP[,c("NAME","R")]
YLD_PERCENTILES_R = YLD_PERCENTILES_R[order(YLD_PERCENTILES_R$M, decreasing = T),]

YLD_PERCENTILES_W = YLD_PERCENTILES_CROP[,c("NAME","W")]
YLD_PERCENTILES_W = YLD_PERCENTILES_W[order(YLD_PERCENTILES_W$W, decreasing = T),]


rm(TONNES_PRES)
rm(t)

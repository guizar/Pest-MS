###
# Fig S4 
# 2 to 4 degree bar for IPM for each crop for each UN region.  for Phi = 0.001.
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
SUM_REGION_IPM = dplyr::select(SUM_REGION_IPM,-fact)
SUM_REGION_IPM = dplyr::filter(SUM_REGION_IPM,phi=="2")
SUM_REGION_IPM = dplyr::select(SUM_REGION_IPM,-phi)

colnames(SUM_REGION_IPM) = c("region","mean","median","min","max","Temp anomaly")

# arrange by LAT
fgLAT = dplyr::select(ALL_2c, one_of("LAT","region"))   # data XY
fgLAT = fgLAT %>%
  group_by(region) %>% 
  summarise_each(funs(mean)) %>% 
  arrange(desc(LAT))

SUM_REGION_IPM$region = as.factor(SUM_REGION_IPM$region)
SUM_REGION_IPM$region = factor(SUM_REGION_IPM$region, levels =rev(fgLAT$region))

# //// PLOT

p = ggplot(SUM_REGION_IPM,aes(x = median, y =region, shape=`Temp anomaly`))
p = p + geom_point()
p = p + geom_line(aes(group=region))
p = p + xlab(label = "Median IPM") + ylab(label = "")

# update theme
p = p %+% mygg
p = p + theme(legend.position="bottom")
p

# Save plot
ppi = 300
plotname = file.path(wdpng,paste("Supp Figure 4",".png",sep = ""))
png(filename=plotname,width=10*ppi, height=8*ppi, res=ppi )
p
dev.off()



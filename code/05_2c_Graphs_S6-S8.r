# /////
#Fig S5 -  S7
# Current yld per HA (YLD) = YLD_HA_X [OK]
# Current YLOSS to insects pests = CLF_X [OK]
# projected additional loss (% of current yield) due to climate change = IYCC_XX [OK]

# GGPLOT OBJECT (MAP) STRUCTURE 
# variables = LON, LAT, phi, crop
# faceting vars = iycc, clf, yld_ha

# --- _2c ---
# Load data and initial variables ------

wd = "~/R/Pest-MS/"
wdpng = "~/R/Pest-MS/png"
wdtables = "~/R/Pest-MS/tables"
wddata = "~/R/Pest-MS/data/"
wdrdata = "~/R/Pest-MS/RData/"
wdfun = "~/R/Pest-MS/fun"

# Load data
load(file.path(wdrdata,"ALL_2c.RData"))

# clean space
source(file.path(wdfun,"clean_space.r")) 


# load libraries used to produce ALL graphs
library(reshape2)
library(classInt)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(R.matlab)
library(arrayhelpers)
library(plyr)
library(tidyr)

# Load functions
source(file.path(wdfun,"cols_gg.r")) 
source(file.path(wdfun,"multiplots_gg.r"))
source(file.path(wdfun,"custom_cut.r"))


# --- PRESENT YIELD PER HA ----
YLD_HA_PRES = ALL_2c[,c("LON","LAT","YLD_HA_M","YLD_HA_R","YLD_HA_W")]
YLD_HA_PRES = melt(YLD_HA_PRES, id.vars=c("LON","LAT"),na.rm=T)

# split PRES cols by CROP
t1 = colsplit(YLD_HA_PRES$variable,"YLD_HA_", c("YLD_HA","crop"))
YLD_HA_PRES$crop = t1$crop
YLD_HA_PRES = YLD_HA_PRES[,-which(colnames(YLD_HA_PRES) %in% "variable")]
rm(t1)


# --- PRESENT YIELD LOST DUE TO PESTS ----
CLF_PRES = ALL_2c[,c("LON","LAT","CLF_M","CLF_R","CLF_W")]
CLF_PRES = melt(CLF_PRES, id.vars=c("LON","LAT"),na.rm=T)

# split PRES cols by CROP
t1 = colsplit(CLF_PRES$variable,"CLF_", c("CLF","crop"))
CLF_PRES$crop = t1$crop
CLF_PRES = CLF_PRES[,-which(colnames(CLF_PRES) %in% "variable")]
rm(t1)


# -- IYCC: Percent of Yield loss do to climate change ----
# Select relevant columns
IYCC_2c = dplyr::select(ALL_2c,LON,LAT,IYCC_M2,IYCC_M3,IYCC_M4,IYCC_R2,IYCC_R3,IYCC_R4,IYCC_W2,IYCC_W3,IYCC_W4)

# Produce long version of the tables
# Change from wide to long
library(reshape2)
IYCC_2c = melt(IYCC_2c, id.vars=c("LON","LAT"),na.rm=T)

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

# -- SUMMARIZE tables ----
IYCC_2c = ddply(IYCC_2c, c("LON","LAT","crop","phi"), summarise,
                iycc= mean(value,na.rm = T))

CLF_PRES = ddply(CLF_PRES, c("LON","LAT"), summarise,
                 clf= mean(value,na.rm = T))

YLD_HA_PRES = ddply(YLD_HA_PRES, c("LON","LAT"), summarise,
                    yld_ha= mean(value,na.rm = T))

# use Phi 3 for IYCC_2c
IYCC_2c = dplyr::filter(IYCC_2c,phi==3)
IYCC_2c = dplyr::select(IYCC_2c,-phi)

# ADD VARIABLE COL = VERY BAD IDEA
# IYCC_2c$variable = "iycc"
# CLF_PRES$variable = "clf"
# YLD_HA_PRES$variable = "yld_ha"

# MERGE
MAP = merge(IYCC_2c,CLF_PRES,by = c("LON","LAT"))
MAP = merge(MAP,YLD_HA_PRES,by = c("LON","LAT"))

# RESHAPE for ggplot
MAP = melt(MAP,id.vars = c("LON","LAT","crop"))

# as factor
MAP$crop = as.factor(MAP$crop)
levels(MAP$crop) = c("Maize", "Rice", "Wheat")

# MAP$phi = as.factor(MAP$phi)
# levels(MAP$phi) = c("Phi 0.01", "Phi 0.001", "Phi 0.0001")


# --- Define distributional range (to remove outliers) ----
wsk_iycc = boxplot.stats(MAP$value[MAP$variable=="iycc"])$stats
wsk_clf = boxplot.stats(MAP$value[MAP$variable=="clf"])$stats
wsk_yld_ha = boxplot.stats(MAP$value[MAP$variable=="yld_ha"])$stats

#subset
MAP_iycc = subset(MAP, variable =="iycc" & value <= wsk_iycc[5] & value >= wsk_iycc[1])
MAP_clf = subset(MAP, variable =="clf" & value <= wsk_clf[5] & value >= wsk_clf[1])
MAP_yld_ha = subset(MAP, variable =="yld_ha" & value <= wsk_yld_ha[5] & value >= wsk_yld_ha[1])

# change percent values
MAP_clf$value = MAP_clf$value *100
MAP_iycc$value = MAP_iycc$value * 100


# --- PLOT ----
# generate basemap
source(file.path(wdfun,"basemap.r"))

# MAP_yld_ha
m1 = ggplot()  + geom_polygon(data=map,aes(x=long, y=lat,group=group), fill="gray", color="gray",size=0.2)
m1 = m1 + geom_raster(data=MAP_yld_ha,aes(fill=value,x=LON, y=LAT),interpolate = T)
m1 = m1 + facet_grid(crop~.)
m1 = m1 + xlab("") + ylab("")
m1 = m1 + scale_fill_gradientn(colours = rev(pal),
                          name="Tonnes/HA")
m1 = m1 %+% mygg
m1 = m1 + theme(legend.position="bottom")
m1

# MAP_clf
m2 = ggplot()  + geom_polygon(data=map,aes(x=long, y=lat,group=group), fill="gray", color="gray",size=0.2)
m2 = m2 + geom_raster(data=MAP_clf,aes(fill=value,x=LON, y=LAT),interpolate = T)
m2 = m2 + facet_grid(crop~.)
m2 = m2 + xlab("") + ylab("")
m2 = m2 + scale_fill_gradientn(colours = rev(pal),
                               name="Yield loss to insect pests (%)")
m2 = m2 %+% mygg
m2 = m2 + theme(legend.position="bottom")
m2

# MAP_iycc
m3 = ggplot()  + geom_polygon(data=map,aes(x=long, y=lat,group=group), fill="gray", color="gray",size=0.2)
m3 = m3 + geom_raster(data=MAP_iycc,aes(fill=value,x=LON, y=LAT),interpolate = T)
m3 = m3 + facet_grid(crop~.)
m3 = m3 + xlab("") + ylab("")
m3 = m3 + scale_fill_gradientn(colours = rev(pal),
                               name="Additional loss (%)\ndue to pest increases")
m3 = m3 %+% mygg
m3 = m3 + theme(legend.position="bottom")
m3

# remove facet names for m1 and m2
m1 = m1 + theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())
m2 = m2 + theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())

# grid plot
grid.arrange(m1,m2,m3, ncol = 3, nrow = 1, widths = c(2,2,2), height = c(1, 1))


# Save plot
ppi = 300
plotname = file.path(wdpng,paste("Supp Figures  S5 -  S7 2c",".png",sep = ""))
png(filename=plotname,width=18*ppi, height=11*ppi, res=ppi )
grid.arrange(m1,m2,m3, ncol = 3, nrow = 1, widths = c(2,2,2), height = c(1, 1))
dev.off()




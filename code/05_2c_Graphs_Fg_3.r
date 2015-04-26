# ///
# Bubble graph: Predicted median increase in insect pest pressure on crops as a function of crop yield (median yield for each country), for A) Maize,  B) Rice, and C) Wheat
# ///

rm(list = ls())

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

# Load functions
source(file.path(wdfun,"cols_gg.r")) 
source(file.path(wdfun,"multiplots_gg.r"))
source(file.path(wdfun,"custom_cut.r"))

# options(scipen=5)

#--- Mean absolute latitude (COLORS) ----
fgLAT = dplyr::select(ALL_2c, one_of("LAT","NAME"))   # data XY
fgLAT = fgLAT %>%
  group_by(NAME) %>% 
  summarise_each(funs(mean)) %>% 
  arrange(desc(LAT))

### add bins of ABS LAT (range)
fgLAT$LAT = abs(fgLAT$LAT) 
fgLAT$LAT = round(fgLAT$LAT) 

## breaks
brk = c(0, 10, 20, 30, 40,max(fgLAT$LAT))

## generate labs
lab = NULL
for (n in seq(length(brk)-1)){
  lb = paste(brk[n]," to ",brk[n+1],sep="")
  lab = c(lab,lb)
}

# add label column
fgLAT = transform(fgLAT, range = cut(LAT,breaks = brk,labels = lab))
fgLAT$range[is.na(fgLAT$range)] = lab[1]
  

#--- TOTAL Y HA (X-axis) ----
fgHA = dplyr::select(ALL_2c, one_of("YLD_HA_M","YLD_HA_R","YLD_HA_W","NAME")) 

# Change from wide to long
fgHA = melt(fgHA, id.vars=c("NAME"),na.rm=T)

# split HA cols by CROP
t1 = colsplit(fgHA$variable,"YLD_HA_", c("YLD_HA","crop"))

# join new columns and delete unused 
fgHA$crop = t1$crop
fgHA = dplyr::select(fgHA,-variable)

## summarize by "NAME","crop"
fgHA = ddply(fgHA, c("NAME","crop"), summarise,
                          yld_ha= median(value,na.rm = T))


#--- YLD TOTAL (TONNES) (circle size)
fgTONNES <- ddply(TONNES_PRES, c("NAME","crop"), summarise,
                                 tonnes= sum(value,na.rm = T))

# max per crop
fgTONMAX <- ddply(fgTONNES, c("crop"), summarise,
                  max= sum(tonnes,na.rm = T))

#--- IPM (Y axis) ----
# //////
# QUESTION: by MEDIAN you are referring to median PHI2 or median across ALL PHIS? below i'm doing MEDIAN ALL
# /////

fgIPM = ddply(DAT_2c,c("fact","phi","NAME","crop"),summarise,
               ipm= median(value,na.rm = T))

fgIPM$ipm = fgIPM$ipm * 100

# remove factors
fgIPM = dplyr::filter(fgIPM,fact =="IPM")
fgIPM = dplyr::select(fgIPM,-fact)

# define distribution range (for outliers)
wsk = boxplot.stats(fgIPM$ipm)$stats
boxplot(fgIPM$ipm)

#--- MERGE DATASETS ----
fg = merge(fgLAT,fgHA,by = "NAME")
fg = merge(fg,fgTONNES,by = c("NAME","crop"))
fg = merge(fg,fgIPM,by = c("NAME","crop"))

# factorize
fg$range = as.factor(fg$range)
fg$crop = as.factor(fg$crop)
levels(fg$crop) = c("Maize", "Rice", "Weat") # rename
fg$phi = as.factor(fg$phi)
levels(fg$phi) = c("Phi 0.01", "Phi 0.001", "Phi 0.0001")

# set breaks/labes for bubble sized
brk = c(3000000,30000000,100000000,240000000)
lab = c("3 Millions","30 Millions","100 Millions","240 Millions")

#--- PLOT ----
p = ggplot(subset(fg, ipm <= wsk[5] & ipm >= wsk[1]), # outliers out
           aes(x = yld_ha, y = ipm)) 
p = p + geom_point(aes(size = tonnes, colour = range),alpha=0.5)

p = p + scale_size_area(breaks=brk, 
                        labels=lab,
                        "Total yield (tones/yr)", 
                        max_size=20)

p = p + scale_color_brewer(palette="Set1", 
                      name="Mean absolute\nLatitude")

p = p + facet_grid(phi~crop)
p = p + guides(colour = guide_legend(override.aes = list(size=8)))
p = p + xlab("Median yield per hectare") + ylab("Increase in pest pressure (%)")
p = p + xlim(c(0,10))
p = p + geom_text(data=subset(fg, tonnes>50000000), aes(label=NAME), size=4)

# update theme
p = p %+% mygg
p

# Save plot
ppi = 300
plotname = file.path(wdpng,paste("Figure 3 2c",".png",sep = ""))
png(filename=plotname,width=14*ppi, height=10*ppi, res=ppi )
p
dev.off()

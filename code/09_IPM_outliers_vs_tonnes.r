# PENDING = percentile rank tonnes

# //// 
# Compare 2c IPM outliers cells against high yielding cells
# This code is continuation of  07_ALLc_Summary_IPM -  need to run it first
# Run the script and look for the following object: SUMMARY_OUTLIERS
# ///

# ///
# SUMMARY_OUTLIERS metadata
# phi = phi associated for each ipm, if NA means that no forecast value was included in this observation
# ipm = ipm value for each crop, if NA means that no forecast value was included in this observation
# tonnes = present yld value in tonnes associated with this observation
# percentile = tonnes percentile rank (each crop has its own rank)
# country_tonnes = present tonnes/yr value associated with this country (for each crop)
# country_percentile =  percentile rank associated with this country (for each crop)
# ////

rm(list = ls())

# Load data and initial variables ------

wd = "~/R/Pest-MS/"
wdpng = "~/R/Pest-MS/png"
wdtables = "~/R/Pest-MS/tables"
wddata = "~/R/Pest-MS/data/"
wdrdata = "~/R/Pest-MS/RData/"
wdfun = "~/R/Pest-MS/fun"

# Load data
load(file.path(wdrdata,"ALL_2c.RData"))
load(file.path(wdrdata,"ALL_4c.RData"))

# clean space
source(file.path(wdfun,"clean_space.r")) 

# libs
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

# functions
source(file.path(wdfun,"cols_gg.r")) 
source(file.path(wdfun,"multiplots_gg.r"))
source(file.path(wdfun,"custom_cut.r"))


# --- 1) Prepare DF (XY coord) with total yield (tonnes) ----
TONNES_XY = ALL_2c[,c("LON","LAT","NAME","region","YLD_TOT_M","YLD_TOT_R","YLD_TOT_W")]

# Change from wide to long
TONNES_XY = melt(TONNES_XY, id.vars=c("LON","LAT","NAME","region"),na.rm=T)

# split PRES cols by CROP
t1 = colsplit(TONNES_XY$variable,"_", c("YIELD_TOT","crop"))
t1 = colsplit(t1$crop,"_", c("YIELD_TOT","crop"))
TONNES_XY$crop = t1$crop
TONNES_XY = TONNES_XY[,-which(colnames(TONNES_XY) %in% "variable")]
rm(t1)

# rename
TONNES_XY$tonnes = TONNES_XY$value
TONNES_XY = TONNES_XY[,-which(colnames(TONNES_XY) %in% "value")]

# add tonnes percentiles
TONNES_XY = TONNES_XY %>%
  group_by(crop) %>%
  mutate(percentile = percent_rank(tonnes))

# --- 2) Compare 2c IPM outliers cells against high yielding cells ----
## 2a) load data
load(file.path(wdrdata,"IPM_2c_XY.RData")) ## IPM spatial points
YLD_CROP_STATS = read.csv( # Total yield and percentiles per country
  file.path(wdtables,"YLD_CROP_STATS_LONG.csv"),header = T)

# rename IPM_2c_XY value column
IPM_2c_XY$ipm = IPM_2c_XY$value
IPM_2c_XY$value = NULL
IPM_2c_XY$fact = NULL

# rename YLD_CROP_STATS vars
YLD_CROP_STATS$country_tonnes = YLD_CROP_STATS$tonnes
YLD_CROP_STATS$country_percentile = YLD_CROP_STATS$percentiles
YLD_CROP_STATS$percentiles = NULL
YLD_CROP_STATS$tonnes = NULL

## 2b) summarize IPM outliers
stats.IPM = boxplot.stats(IPM_2c_XY$ipm)

# extract outliers
IPM_2c_OUTLIERS = IPM_2c_XY[IPM_2c_XY$ipm  %in% stats.IPM$out[stats.IPM$out>1],]

## 2c) subset outliers from TONNES object
SUMMARY_OUTLIERS = dplyr::semi_join(TONNES_XY, IPM_2c_OUTLIERS, by = c("LON","LAT"))

## 2d) MERGE DATASETS
# join SUMMARY_OUTLIERS and IPM_2c_OUTLIERS data (add ipm column)
SUMMARY_OUTLIERS = dplyr::right_join(IPM_2c_OUTLIERS, SUMMARY_OUTLIERS, by = c("LON","LAT","NAME","region","crop"))

# join SUMMARY_OUTLIERS with YLD_CROP_STATS data (add tonnes as  country_tonnes, and percentiles as country_percentile)
SUMMARY_OUTLIERS = dplyr::right_join(YLD_CROP_STATS, SUMMARY_OUTLIERS, by = c("NAME","crop"))

# horizontal re-ordering (columns)
cols = c("LON","LAT","NAME","region","crop","phi","ipm","tonnes","percentile","country_tonnes","country_percentile")
SUMMARY_OUTLIERS = SUMMARY_OUTLIERS[,cols]
rm(cols)

# write.csv
write.csv(SUMMARY_OUTLIERS, file=file.path(wdtables,"SUMMARY_IPM_OUTLIERS.csv"),row.names=F)


# --- 3) PLOT SUMMARY_OUTLIERS data ----
gg = SUMMARY_OUTLIERS
options(scipen=99)

# IPM outliers (IPM > 1), tonnes and IPM values per grid cell
# x = ipm
# y = countries arranged by total tonnes/yr
# colors = tonnes percentile rank (each crop has its own rank)
# each pannel represent a different phi

gg = dplyr::filter(gg,phi != "NA")

fg = ggplot(data = gg, 
            aes(x=ipm, 
                y=reorder(NAME, country_tonnes), 
                color=percentile, shape=crop), 
            group=NAME)
fg = fg + geom_point(position = position_identity(),alpha=0.7)
fg = fg + scale_color_gradientn(colours = rev(pal))
fg = fg + facet_grid(phi~.)
fg = fg + ggtitle("IPM outliers (IPM > 1)\ntonnes and IPM values per grid cell")
fg = fg + ylab("")
fg
ggsave(filename = file.path(wdpng,"ipm_outliers_scale_IPM_full.png"),scale = 1.5)

fg = fg + xlim(c(0,2000))
fg
ggsave(filename = file.path(wdpng,"ipm_outliers_scale_IPM_2K.png"),scale = 1.5)

## Density plots for percentile ranks
fg = ggplot(gg, aes(x=percentile)) 
fg = fg + geom_histogram(aes(y=(..count../sum(..count..))*100,fill=crop,alpha=0.7))
fg = fg + geom_density(data = gg, aes(x=percentile))
fg = fg + ggtitle("IPM outliers\n Tonnes percentile rank by grid cell")
fg = fg + ylab("IPM count (re-scaled)")
fg
ggsave(filename = file.path(wdpng,"ipm_percentile_ranks.png"),scale = 1.5)


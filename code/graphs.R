# Load data and initial variables

wd = "~/R/Pest-MS/"
wdpng = "~/R/Pest-MS/png"
wddata = "~/R/Pest-MS/data"
wdfun = "~/R/Pest-MS/functions"

# # functions
# source(file.path(wdfun,"color bar.R"))
# source(file.path(wdfun,"summarySE function.R"))
# source(file.path(wdfun,"multiplot.R"))


# Load ALL.CSV --- I assume this is 2c
setwd(wddata)
ALL<-read.csv("ALL.csv", header=TRUE,na.strings = c("#VALUE!", "#N/A", "N/A", "NA", ""))



# //////////////////////// FIGURE 1) 

# Fig 1a) Map of Fractional change in population metablism

# which variables I'm using?
# METABOLISM GRAPHIC = MET_AVG3
# POPULATION GRAPHIC = POP_AVG3
# IPM AVERAGE ACROSS ALL THREE CROPS = IPM_AVG3

# subset data  
fg = ALL[,c("LAT","LON","MET_AVG3","POP_AVG3","IPM_AVG3")]

# change format from wide to long (for ggplot)
library(reshape2)
fg = melt(fg, id.vars=c("LAT","LON"))

## Generate breaks
library(classInt)
breaks= c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)
fg$brks <- cut(fg$value, breaks =breaks, labels = as.character(breaks[-1]), right = F)

## change "variable" names [for ggplot]
levels(fg$variable) = c("Metabolism", "Population", "Total")

## Prepare color palette
library(RColorBrewer)
pal <- brewer.pal("Spectral",n = 11)
cols<- colorRampPalette(pal)

# generate basemap
library(ggplot2)
library(mapdata)
map = map_data("world")

# out antartica and big lakes
map = map[!map$region=="Antarctica",]
lakes1 = grep(unique(map$region),pattern = "Lake",value = T,ignore.case = T)
lakes2 = grep(unique(map$subregion),pattern = "lake",value = T,ignore.case = T)
for (i in lakes1){
  map = map[!map$region==i,]	
}

rm(lakes1)
rm(lakes2)
rm(i)

## plot
m1 = ggplot()  + geom_polygon(data=map,aes(x=long, y=lat,group=group), fill="gray", color="gray",size=0.2)
m1 = m1 + geom_raster(data=fg,aes(fill=brks,x=LON, y=LAT),interpolate = T)
m1 = m1 + scale_fill_manual(values =rev(cols(12)))
m1 = m1 +  ggtitle("Fractional change in population metabolism")
m1 = m1 + facet_wrap(~variable, ncol = 1,  scales="free")
m1 = m1 + guides(fill=guide_legend(title=NULL))

# add theme to map
m1 = m1 + theme_bw()
m1 = m1 + xlab("") + ylab("")

m1

# Save plot
# ppi = 300
# plotname = file.path(wdpng,paste("fig1",".png",sep = ""))
# png(filename=plotname,width=9*ppi, height=15*ppi, res=ppi )
# m1
# dev.off()

# Fig 1b) 3 adjacent panels indicating phi values per UN region
fg1b = ALL[,c("region","MET_M2","MET_R2","MET_W2","MET_M3","MET_R3","MET_W3","MET_M4","MET_R4","MET_W4","POP_M2","POP_R2","POP_W2","POP_M3","POP_R3","POP_W3","POP_M4","POP_R4","POP_W4","IPM_M2","IPM_R2","IPM_W2","IPM_M3","IPM_R3","IPM_W3","IPM_M4","IPM_R4","IPM_W4")]


# a) create a variable column. this will be used to apply the colsplit() fun to separate by crops, phi and factor

# Change from wide to long
library(reshape2)
fg1b = melt(fg1b, id.vars="region")

# split cols
t = colsplit(fg1b$variable,"_", c("fact","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi"))
t = colsplit(fg1b$variable,"_", c("fact","crop"))
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi"))

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

fg1b = fg1b[,-which(names(fg1b) %in% "variable")]
fg1b$crop = t$crop
fg1b$phi = t$phi
fg1b$fact = t$fact

rm(t)
rm(c2)
rm(c1)

# b) Summarize and prepare data for ggplot:

# sumarize
library(plyr)
fg1b = ddply(fg1b, c("region","crop","phi","fact"), summarise,value  = mean(value, na.rm=T))

# change colnames
fg1b = dcast(fg1b, region + crop + fact ~ phi, value.var="value")
names(fg1b)[names(fg1b)=="2"] <- "phi2"
names(fg1b)[names(fg1b)=="3"] <- "phi3"
names(fg1b)[names(fg1b)=="4"] <- "phi4"

# now melt phi 
fg1b = melt(fg1b, id.vars=c("region","crop","fact"))
fg1b = fg1b[!is.nan(t$value),]

# prepare for plotting
# as.fcators (sets facetting order)
fg1b$fact = factor(fg1b$fact, levels=c('MET','POP','IPM'))

# change crop names
fg1b$crop = as.factor(fg1b$crop )
levels(fg1b$crop) = c("Maize", "Rice", "Wheat")

# PLOT
library(ggplot2)
p = ggplot(fg1b, aes(x = region, y = value)) + geom_boxplot() 
p = p + coord_flip()
p = p + facet_grid(fact~crop) 
p = p + ylim(ymin=0, ymax=0.5)
p = p + xlab(label = "") + ylab(label = "Phi range") 
p

# Save plot
ppi = 300
plotname = file.path(wdpng,paste("fig1b",".png",sep = ""))
png(filename=plotname,width=10*ppi, height=10*ppi, res=ppi )
p
dev.off()

## MISSING: ALLGIN fg1 (maps) and fg1b

# //////////////////////// 


# //////////////////////// 


###------- FIGURE S1:  Maps of demographic change for each crop, for each value of Phi (9 maps) showing consistency of pattern across models.  

# which variables I'm using?
# IPM_xy
  
fg = ALL[,c("LAT","LON","IPM_M2","IPM_M3","IPM_M4","IPM_R2","IPM_R3","IPM_R4","IPM_W2","IPM_W3","IPM_W4")]




  
### --------- FIGURE S2 : latitudinal median change in insect pest pressure at high (black), medium, (yellow) and low (blue) overwinter survival for maize (A), Rice (B), and Wheat (C).

# which variables am I using?
# IMP_xy


fg = ALL[,c("LAT","LON","IPM_M2","IPM_M3","IPM_M4","IPM_R2","IPM_R3","IPM_R4","IPM_W2","IPM_W3","IPM_W4")]



## ---------- FIGURE 3: Predicted median increase in insect pest pressure on crops as a function of crop yield (median yield for each country), for A) Maize,  B) Rice, and C) Wheat BUBBLE PLOTS

setwd(wddata)

# read csvs
MAIZE_c<-read.csv("MAIZE by country.csv", header=TRUE)
MAIZE_c<-MAIZE_c[complete.cases(MAIZE_c),]
RICE_c<-read.csv("RICE by country.csv", header=TRUE)
RICE_c<-RICE_c[complete.cases(RICE_c),]
WHEAT_c<-read.csv("WHEAT by country.csv", header=TRUE)
WHEAT_c<-WHEAT_c[complete.cases(WHEAT_c),]

# subset = countries with total yield = .001 of top yield
MAIZE_c<-subset(MAIZE_c,TOT_YLD_TOT_M >= max(TOT_YLD_TOT_M*.001)) 
RICE_c<-subset(RICE_c,TOT_YLD_TOT_M >= max(TOT_YLD_TOT_M*.001)) 
WHEAT_c<-subset(WHEAT_c,TOT_YLD_TOT_M >= max(TOT_YLD_TOT_M*.001)) 

# what variables do i NEED?
# TOT_YLD_TOT_x
# WMEAN_YLD_HA_x
# WMEAN_IPM_xy
# MEAN_LAT_x


####
# MAIZE_c = MAIZE_c$TOT_YLD_TOT_M
# 
# 
# ggplot(asd_data, aes(x=prev2003, y=asd_diff, weight=denom2003, colour=octile, size=denom2003)) + 
#   geom_point( alpha=0.8, guide="none") + 
#   scale_size_area(breaks=c(250, 500, 1000, 10000, 50000), "2002 District\nElementary School\nPopulation", max_size=20) + 
#   stat_smooth(method="rlm", size=0.5, colour="black", alpha=0.4, level=0.95)+
#   scale_colour_brewer(palette="Spectral", type="qual",name="2002 Autism\nPrevalence Octile") + 
#   coord_equal(ratio=1/2)+
#   guides(colour = guide_legend(override.aes = list(alpha = 1)))+
#   ggtitle("Figure 4. Change in Autism Prevalence between 2002 and 2008 vs Baseline (2002) Prevalence,\n 
#        Wisconsin Elementary School Districts (with weighted linear best-fit line and 95% confidence band)") +
#   scale_x_continuous("2002 Autism Prevalence (per 1,000)") + 
#   scale_y_continuous("Change in Autism Prevalence (per 1,000) between 2002 and 2008")
####










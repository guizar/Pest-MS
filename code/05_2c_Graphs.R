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



# Fig 1a) Map of Fractional change in population metablism ----
# Variables
# METABOLISM GRAPHIC = MET_AVG3
# POPULATION GRAPHIC = POP_AVG3
# IPM AVERAGE ACROSS ALL_2c THREE CROPS = IPM_AVG3

# subset data  
fg = ALL_2c[,c("LAT","LON","MET_AVG3","POP_AVG3","IPM_AVG3")]

# change format from wide to long (for ggplot)
fg = melt(fg, id.vars=c("LAT","LON"))

## Generate breaks
breaks= c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)
fg$brks <- cut(fg$value, breaks =breaks, labels = as.character(breaks[-1]), right = F)

## change "variable" names [for ggplot]
levels(fg$variable) = c("Metabolism", "Population", "Total")

# generate basemap
source(file.path(wdfun,"basemap.r"))

## plot
m1 = ggplot()  + geom_polygon(data=map,aes(x=long, y=lat,group=group), fill="gray", color="gray",size=0.2)
m1 = m1 + geom_raster(data=fg,aes(fill=brks,x=LON, y=LAT),interpolate = T)
m1 = m1 + scale_fill_manual(values =rev(cols(12)))
m1 = m1 + theme(title=element_blank())
m1 = m1 + facet_grid(variable~.)
m1 = m1 + guides(fill=guide_legend(title=NULL))
m1 = m1 + theme(panel.grid = element_blank())
m1 = m1 + xlab("") + ylab("")

# update theme
m1 = m1 %+% mygg
m1

# Save plot
ppi = 300
plotname = file.path(wdpng,paste("Figue 1A 2c",".png",sep = ""))
png(filename=plotname,width=9*ppi, height=15*ppi, res=ppi )
m1
dev.off()

# Fig 1b) Summary of fractional change (horizontal plot, phi3 in the middle) ----
# Variables:
# MET | POP | IPM
# regions
# crops
# phi
# Produce this graph so that hey are arranged by latitudinal central tendencies (mean latitude)

# data
fg1b = DAT_2c # data
fgXY = dplyr::select(ALL_2c, one_of("LON","LAT","region"))   # data XY

# get "mean latitude" and arrnage by decreasing LATs
fgXY= fgXY %>%
  group_by(region) %>% 
  summarise_each(funs(mean)) %>% 
  arrange(desc(LAT))

# sumarize data
fg1b = dplyr::select(fg1b,-NAME) # remove country names for this summary

fg1b = fg1b %>%
  group_by(region,crop,phi,fact) %>% 
  summarise_each(funs(mean(.,na.rm=T)))

# re-arrange data: add columns for each phi
fg1b = dcast(fg1b, region + crop + fact ~ phi, value.var="value")
names(fg1b)[names(fg1b)=="2"] <- "phi2"
names(fg1b)[names(fg1b)=="3"] <- "phi3"
names(fg1b)[names(fg1b)=="4"] <- "phi4"

# set facetting order
fg1b$fact = factor(fg1b$fact, levels=c('MET','POP','IPM'))

# change crop names
fg1b$crop = as.factor(fg1b$crop )
levels(fg1b$crop) = c("Maize", "Rice", "Wheat")

# change vertical arrangement by mean lat
fg1b$region = as.factor(fg1b$region)
levels(fg1b$region) = rev(fgXY$region)

# //// PLOT

# init variables
pd = position_dodge(width = 1)

p = ggplot(fg1b, aes(x = region, y =phi3, color=crop)) 
p = p + geom_errorbar(aes(ymin=phi2, ymax=phi4), width=0.05, position=pd, width = 0.5,stat="identity")
p = p + geom_point(position = pd)
p = p + coord_flip()
p = p + facet_grid(fact~.) 
p = p + ylim(xmin=-0.25, ymax=0.5)
p = p + xlab(label = "") + ylab(label = "Fractional change")

# horizontal lines
p = p + geom_vline(xintercept=seq(0.5, length(unique(fg1b$region)), 1), lwd=0.2, colour="black")
p = p + theme(panel.grid.major.y = element_blank())
p = p + theme(title=element_blank())

# scales
labs = levels(fg1b$crop)
p = p + scale_color_manual(name = "Crops", labels = labs, breaks=labs, values=palcrop)

# update theme
p = p %+% mygg
p

# # Save plot
ppi = 300
plotname = file.path(wdpng,paste("Figue 1B 2c",".png",sep = ""))
png(filename=plotname,width=6*ppi, height=8*ppi, res=ppi )
p
dev.off()
# ggsave(plotname,p,scale = 2,dpi=150)


# Fig 1a) and 1b) combined ----

# plot m1 again (no ggtitle)
m1 = ggplot()  + geom_polygon(data=map,aes(x=long, y=lat,group=group), fill="gray", color="gray",size=0.2)
m1 = m1 + geom_raster(data=fg,aes(fill=brks,x=LON, y=LAT),interpolate = T)
m1 = m1 + scale_fill_manual(values =rev(cols(12)))
m1 = m1 + facet_grid(variable~.)
m1 = m1 + guides(fill=guide_legend(title=NULL))
m1 = m1 + theme(panel.grid = element_blank())
m1 = m1 %+% mygg
m1 = m1 + theme(title=element_blank())

# remove unused stuff m1
m1 = m1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank())

# remove xlab p
p = p + theme(axis.title.x = element_blank())

# place legends at the bottom (1a and 1b)
m1 = m1 + theme(legend.position="bottom")
p = p + theme(legend.position="bottom")

# remove facet names for p
p = p + theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())

# grid plot
grid.arrange(p,m1, ncol = 2, nrow = 1, widths = c(1,2), height = c(1, 1))
# main = textGrob("Fractional change in population metabolism",gp=gpar(fontsize=20,font=2))

# Save plot
ppi = 300
plotname = file.path(wdpng,paste("Figure 1 2c",".png",sep = ""))
png(filename=plotname,width=11*ppi, height=10*ppi, res=ppi )
grid.arrange(p,m1, ncol = 2, nrow = 1, widths = c(1,1.5), height = c(1, 1))
dev.off()





# --- Fig 2 (smooth): Globally integrated loss of crop yield  ----
# x2 = temp anomaly
# y2 = yield loss
#
# [1:4] = GFDL, MPI, IPSL, Hadley
# [1:3] = Wheat, Rice, Maize
# [1:2] = Ph 0.01, Ph 0.0001
# [1:100] =  Temporal anomaly
###

data = readMat(file.path(wddata,"MS_Fig2.mat"))

dfy = array2df(data$y)
dfx = array2df(data$x)
dfx = dfx[c("d1","d2","d3","data$x")]
dfXy = cbind(dfx,dfy$"data$y") 

colnames(dfXy) = c("Model","Crop","Phi","x","y") 

# Change values
# MODELS
dfXy$Model[dfXy$Model==1] ="GFDL"
dfXy$Model[dfXy$Model==2] ="MPI"
dfXy$Model[dfXy$Model==3] ="IPSL"
dfXy$Model[dfXy$Model==4] ="Hadley"

# CROPS
dfXy$Crop[dfXy$Crop==1] ="Wheat"
dfXy$Crop[dfXy$Crop==2] ="Rice"
dfXy$Crop[dfXy$Crop==3] ="Maize"

# phi
dfXy$Phi[dfXy$Phi==1] ="Phi 0.0001"
dfXy$Phi[dfXy$Phi==2] ="Phi 0.01"

# Rescale Y
dfXy$y =  dfXy$y/1000000

# change crop names
dfXy$Crop = as.factor(dfXy$Crop )
levels(dfXy$Crop) = c("Maize", "Rice", "Wheat")

library (ggplot2)

# PLOT
p = ggplot(dfXy, aes(x=x, y=y, color=Crop)) 
# p = ggplot(dfXy, aes(x=x, y=y, color=Crop, shape=Phi)) # (plotly)
p = p + geom_smooth(method=lm, size= 1) 
# p = p + geom_point() #  (plotly)
p = p + facet_wrap(~Phi, ncol = 1) # (remove for plotly)
# p + geom_ribbon(alpha=0.2)
# annotations
ylab=expression(bold(Yield~loss~~"(Ton/yr)"*10^"6"))
# ylab= "Yield loss (Ton/yr) *10^6" # (plotly)
p = p + theme(title=element_blank())
p = p + xlab("Temp anomaly") + ylab(ylab)

## Fixed Y 
# p + scale_y_discrete(breaks=seq(0, max(dfXy$y), 3))

# scales
labs = levels(dfXy$Crop)
p = p + scale_color_manual(name = "Crops", labels = labs, breaks=labs, values=palcrop)

# theme
p = p %+% mygg
p = p + theme(legend.position="bottom")
p

# Save plot
ppi = 300
plotname = file.path(wdpng,paste("Figure 2 v1",".png",sep = ""))
png(filename=plotname,width=6*ppi, height=8*ppi, res=ppi )
p
dev.off()


# --- Fig 2 (dotplot): Globally integrated loss of crop yield  ----

data = readMat(file.path(wddata,"MS_Fig2.mat"))

dfy = array2df(data$y)
dfx = array2df(data$x)
dfx = dfx[c("d1","d2","d3","data$x")]
dfXy = cbind(dfx,dfy$"data$y") 

colnames(dfXy) = c("Model","Crop","Phi","x","y") 

# Change values
# MODELS
dfXy$Model[dfXy$Model==1] ="GFDL"
dfXy$Model[dfXy$Model==2] ="MPI"
dfXy$Model[dfXy$Model==3] ="IPSL"
dfXy$Model[dfXy$Model==4] ="Hadley"

# CROPS
dfXy$Crop[dfXy$Crop==1] ="Wheat"
dfXy$Crop[dfXy$Crop==2] ="Rice"
dfXy$Crop[dfXy$Crop==3] ="Maize"

# PH
dfXy$Phi[dfXy$Phi==1] ="Phi 0.0001"
dfXy$Phi[dfXy$Phi==2] ="Phi 0.01"

# Rescale Y
dfXy$y =  dfXy$y/1000000

# change crop names
dfXy$Crop = as.factor(dfXy$Crop )
levels(dfXy$Crop) = c("Maize", "Rice", "Wheat")

# PLOT
p = ggplot(dfXy, aes(x=x, y=y, color=Crop, shape=Phi)) # (plotly)
p = p + geom_jitter(alpha=0.5) #  (for plotly use geom_point)
# no faceting remove for plotly

# annotations
ylab=expression(bold(Yield~loss~~"(Ton/yr)"*10^"6"))
# ylab= "Yield loss (Ton/yr) *10^6" # (plotly)
# p = p + ggtitle("Globally integrated loss of crop yield") 
p = p + theme(title = element_blank())
p = p + xlab("Temp anomaly") + ylab(ylab)

# scales
labs = levels(dfXy$Crop)
p = p + scale_color_manual(name = "Crops", labels = labs, breaks=labs, values=palcrop)

# theme
p = p %+% mygg
p = p + theme(legend.position="bottom")
p

# Save plot
ppi = 300
plotname = file.path(wdpng,paste("Figure 2 v2",".png",sep = ""))
png(filename=plotname,width=8*ppi, height=8*ppi, res=ppi )
p
dev.off()



# --- Fig S1: Maps of demographic change for each crop, for each value of Ph ----
# |IPM|
# Phi
# Crops

fg = dplyr::select(ALL_2c,LAT,LON,IPM_M2,IPM_M3,IPM_M4,IPM_R2,IPM_R3,IPM_R4,IPM_W2,IPM_W3,IPM_W4) 

fg = melt(fg, id.vars=c("LON","LAT"), na.rm=T) 

#### split fg cols by CROPS and PHI
t = colsplit(fg$variable,"_", c("IPM","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## gives you the crop name initials
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## gives you the phi values

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

fg = fg[,-which(names(fg) %in% "variable")]
fg$crop = t$crop
fg$phi = t$phi
rm(t)
rm(c2)
rm(c1)
####

### Generate breaks
# explore breaks
wsk = boxplot.stats(fg$value)$stats
by = 0.05
seq.int(wsk[1]-by,wsk[5],by=by) # play with by 

# generate breaks (see custom_cut.r fir instructions)
fg$brks = custom.cut(fg$value,by = 0.05,digits = 3)

### Change labels (for ggplot)
fg$crop = as.factor(fg$crop)
levels(fg$crop) = c("Maize", "Rice", "Wheat") # rename

fg$phi = as.factor(fg$phi)
levels(fg$phi) = c("Phi 0.01", "Phi 0.001", "Phi 0.0001") # rename


# PLOT --- With breaks 
# generate basemap
source(file.path(wdfun,"basemap.r"))

## plot
m1 = ggplot()  + geom_polygon(data=map,aes(x=long, y=lat,group=group), fill="gray", color="gray",size=0.2) #basemap

m1 = m1 + geom_raster(data=fg,aes(fill=brks,x=LON, y=LAT),interpolate = T)
m1 = m1 +  theme(legend.title=element_blank())
m1 = m1 + facet_grid(crop~phi)
m1 = m1 + guides(fill=guide_legend(title=NULL))
m1 = m1 + theme(panel.grid = element_blank())
m1 = m1 + xlab("") + ylab("")

m1 = m1 + scale_fill_manual(values =cols(length(unique(fg$brks))))

# update theme
m1 = m1 %+% mygg
m1

# Save plot
ppi = 300
plotname = file.path(wdpng,paste("Supp Figure 1 2c",".png",sep = ""))
png(filename=plotname,width=14*ppi, height=8*ppi, res=ppi )
m1
dev.off()
# ---


# --- Fig S2: Latitudinal median change in insect pest pressure at high (black), medium, (yellow) and low (blue) overwinter survival for maize (A), Rice (B), and Wheat (C) ----

# data 
fg = dplyr::select(ALL_2c, one_of("LAT",grep(pattern = "^IPM_[A-Z][0-9]",colnames(ALL_2c),value = T)))

#### split fg cols by CROPS and PHI
fg = melt(fg,id.vars = c("LAT"),na.rm = T)
t = colsplit(fg$variable,"_", c("IPM","crop"))
c2 = colsplit(t$crop,"[0-9]", c("crop","phi")) ## crop 
c1 = colsplit(t$crop,"[A-Z]", c("crop","phi")) ## phi

# join new columns and delete unused 
c1$crop = c2$crop
t$crop = c1$crop
t$phi = c1$phi

fg = fg[,-which(names(fg) %in% "variable")]
fg$crop = t$crop
fg$phi = t$phi
rm(t)
rm(c2)
rm(c1)
####

# summarize medians by LON,LAT

fg = fg %>%
  group_by(LAT,crop,phi) %>% 
  summarise_each(funs(median(.,na.rm=T))) %>%
  arrange(.,LAT)

# fg = fg[order(fg$LAT,decreasing = T),]

## Change labels (for ggplot)
fg$crop = as.factor(fg$crop)
levels(fg$crop) = c("Maize", "Rice", "Wheat") # rename

fg$phi = as.factor(fg$phi)
levels(fg$phi) = c("Phi 0.01", "Phi 0.001", "Phi 0.0001") # rename

# calculate stats for outliers
wsk = boxplot.stats(fg$value)$stats

# //// PLOT
p = ggplot(subset(fg, value <= wsk[5] & value >= wsk[1]), # outliers out
           aes(x = value, y =LAT, color=phi)) 

p = p + facet_grid(.~crop) 

p = p + geom_path (aes(group = phi))

p = p + scale_color_manual(values =palcrop)
p = p + xlab(label = "IPM") + ylab(label = "")
p = p + theme(legend.title=element_blank())

# update theme
# p = p %+% mygg
p

# Save plot
ppi = 300
plotname = file.path(wdpng,paste("Supp Figure 2 2c",".png",sep = ""))
png(filename=plotname,width=10*ppi, height=8*ppi, res=ppi )
p
dev.off()

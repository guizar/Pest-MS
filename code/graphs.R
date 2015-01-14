# Load data and initial variables ------

wd = "~/R/Pest-MS/"
wdpng = "~/R/Pest-MS/png"
wddata = "~/R/Pest-MS/data"
wdfun = "~/R/Pest-MS/functions"

# Load data
setwd(wddata)
ALL = read.csv("ALL.csv", header=TRUE,na.strings = c("#VALUE!", "#N/A", "N/A", "NA", ""))
dat = read.csv(file.path(wddata,"dat.csv")) # data produced by the summarize_by_region.r script 
dat = dat[,-1]

# Load functions
source(file.path(wdfun,"cols_gg.r")) 
source(file.path(wdfun,"multiplots_gg.r"))


# Fig 1a) Map of Fractional change in population metablism ----
# Variables
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

# generate basemap
source(file.path(wdfun,"basemap.r"))

## plot
library(ggplot2)
m1 = ggplot()  + geom_polygon(data=map,aes(x=long, y=lat,group=group), fill="gray", color="gray",size=0.2)
m1 = m1 + geom_raster(data=fg,aes(fill=brks,x=LON, y=LAT),interpolate = T)
m1 = m1 + scale_fill_manual(values =rev(cols(12)))
m1 = m1 +  ggtitle("Fractional change in population metabolism")
m1 = m1 + facet_grid(variable~.)
m1 = m1 + guides(fill=guide_legend(title=NULL))
m1 = m1 + theme(panel.grid = element_blank())
m1 = m1 + xlab("") + ylab("")

# update theme
m1 = m1 %+% mygg
m1

# Save plot
# ppi = 300
# plotname = file.path(wdpng,paste("fig1",".png",sep = ""))
# png(filename=plotname,width=9*ppi, height=15*ppi, res=ppi )
# m1
# dev.off()




# Fig 1b) Summary of fractional change (horizontal plot, phi3 in the middle) ----
# Variables:
# MET | POP | IPM
# regions
# crops
# phi

# data
fg1b = dat
levels(fg1b$region)[1] = "Australia and NZ"

# sumarize
library(plyr)
fg1b = ddply(fg1b, c("region","crop","phi","fact"), summarise,value  = mean(value, na.rm=T)) # mean value / plot base

# re-arrange data: add columns for each phi
library(reshape2)
fg1b = dcast(fg1b, region + crop + fact ~ phi, value.var="value")
names(fg1b)[names(fg1b)=="2"] <- "phi2"
names(fg1b)[names(fg1b)=="3"] <- "phi3"
names(fg1b)[names(fg1b)=="4"] <- "phi4"

# set facetting order
fg1b$fact = factor(fg1b$fact, levels=c('MET','POP','IPM'))

# change crop names
fg1b$crop = as.factor(fg1b$crop )
levels(fg1b$crop) = c("Maize", "Rice", "Wheat")

# //// PLOT
library(ggplot2)

# init variables
pd = position_dodge(width = 1)

p = ggplot(fg1b, aes(x = region, y =phi3, color=crop)) 
p = p + geom_errorbar(aes(ymin=phi2, ymax=phi4), width=0.05, position=pd, width = 0.5) 
p = p + geom_point(position = pd)
p = p + coord_flip()
p = p + facet_grid(fact~.) 
p = p + ylim(xmin=-0.25, ymax=0.5)
p = p + xlab(label = "") + ylab(label = "Fractional change")

# horizontal lines
p = p + geom_vline(xintercept=seq(0.5, length(unique(fg1b$region)), 1), lwd=0.2, colour="black")
p = p + theme(panel.grid.major.y = element_blank())

# scales
labs = levels(fg1b$crop)
p = p + scale_color_manual(name = "Crops", labels = labs, breaks=labs, values=palcrop)

# update theme
p = p %+% mygg
p

# # Save plot
# ppi = 300
# plotname = file.path(wdpng,paste("fig1b",".png",sep = ""))
# png(filename=plotname,width=6*ppi, height=8*ppi, res=ppi )
# p
# dev.off()




# Fig 1c) Latitudinal median change in insect pest pressure ----
 



# Fig 1) 1a and 1b combined ----

# plot m1 again (no ggtitle)
m1 = ggplot()  + geom_polygon(data=map,aes(x=long, y=lat,group=group), fill="gray", color="gray",size=0.2)
m1 = m1 + geom_raster(data=fg,aes(fill=brks,x=LON, y=LAT),interpolate = T)
m1 = m1 + scale_fill_manual(values =rev(cols(12)))
m1 = m1 + facet_grid(variable~.)
m1 = m1 + guides(fill=guide_legend(title=NULL))
m1 = m1 + theme(panel.grid = element_blank())
m1 = m1 %+% mygg

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
library(gridExtra)
# grid.arrange(p,m1, ncol = 2, nrow = 1, widths = c(1,2), height = c(1, 1),main = textGrob("Fractional change in population metabolism",gp=gpar(fontsize=20,font=2)))

# Save plot
# ppi = 300
# plotname = file.path(wdpng,paste("Figure_1",".png",sep = ""))
# png(filename=plotname,width=11*ppi, height=11*ppi, res=ppi )
# grid.arrange(p,m1, ncol = 2, nrow = 1, widths = c(1,1.5), height = c(1, 1),main = textGrob("Fractional change in population metabolism",gp=gpar(fontsize=20,font=2)))
# dev.off()




# Fig 2a) Globally integrated loss of crop yield (smooth) ----
library(R.matlab)
library(arrayhelpers)

setwd(wddata)
data = readMat("MS_Fig2.mat")

#### Data composition
# x2 = temp anomaly
# y2 = yield loss
#
# [1:4] = GFDL, MPI, IPSL, Hadley
# [1:3] = Wheat, Rice, Maize
# [1:2] = Ph 0.01, Ph 0.0001
# [1:100] =  Temporal anomaly
###

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
dfXy$Phi[dfXy$Phi==1] ="Phi 2"
dfXy$Phi[dfXy$Phi==2] ="Phi 4"

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
p = p + ggtitle("Globally integrated loss of crop yield") + xlab("Temp anomaly") + ylab(ylab)

## Fixed Y 
# p + scale_y_discrete(breaks=seq(0, max(dfXy$y), 3))

# scales
labs = levels(dfXy$Crop)
p = p + scale_color_manual(name = "Crops", labels = labs, breaks=labs, values=palcrop)

# theme
p = p %+% mygg
p

# Save plot
# ppi = 300
# plotname = file.path(wdpng,paste("fig2",".png",sep = ""))
# png(filename=plotname,width=7*ppi, height=7*ppi, res=ppi )
# p
# dev.off()


# Fig 2b) Globally integrated loss of crop yield (dotplot) ----

# PLOT

p = ggplot(dfXy, aes(x=x, y=y, color=Crop, shape=Phi)) # (plotly)
p = p + geom_jitter(alpha=0.5) #  (for plotly use geom_point)
# no faceting remove for plotly

# annotations
ylab=expression(bold(Yield~loss~~"(Ton/yr)"*10^"6"))
# ylab= "Yield loss (Ton/yr) *10^6" # (plotly)
p = p + ggtitle("Globally integrated loss of crop yield") + xlab("Temp anomaly") + ylab(ylab)

# scales
labs = levels(dfXy$Crop)
p = p + scale_color_manual(name = "Crops", labels = labs, breaks=labs, values=palcrop)

# theme
p = p %+% mygg
p

# Save plot
# ppi = 300
# plotname = file.path(wdpng,paste("fig2b",".png",sep = ""))
# png(filename=plotname,width=7*ppi, height=7*ppi, res=ppi )
# p
# dev.off()



# Fig 3) 2 to 4 degree bar for IPM for each crop for each UN region.----
# Variables
# | IPM 2 deg| IPM 4 deg | 
# regions
# crops
# Phi = 0.001
# pending: need to get 4 deg values




# Fig S1a) Maps of demographic change for each crop, for each value of Ph####
# |IPM|
# Phi
# Crops

fg = ALL[,c("LAT","LON","IPM_M2","IPM_M3","IPM_M4","IPM_R2","IPM_R3","IPM_R4","IPM_W2","IPM_W3","IPM_W4")]


# Fig S1b) Adjacent panels summarizing fractional change by regions and crops (all data)----
fg = dat

# leave IPM only
fg = fg[fg$fact == "IPM",]

# remove unused cols
fg = fg[,c(-1,-5)]

# change crop names
fg$crop = as.factor(fg$crop )
levels(fg$crop) = c("Maize", "Rice", "Wheat")

# PLOT # need to see how to combine flip and facet_wrap
library(ggplot2)
p = ggplot(fg, aes(x = region, y = value)) + geom_boxplot() 
# p = p + coord_flip()
p = p + facet_grid(.~crop)
p = p + ylim(ymin=-0.5, ymax=0.5)
p = p + xlab(label = "") + ylab(label = "Fractional change") 
p

# Save plot
# ppi = 300
# plotname = file.path(wdpng,paste("fig1b_v2",".png",sep = ""))
# png(filename=plotname,width=10*ppi, height=10*ppi, res=ppi )
# p
# dev.off()



# Fig S3) Relationship between yield (relative yield, in 10 percentile bins) and median change in insect pest pressure ----
# Variables
# IPM
# Phi
# relative yield in 10 percentile bins




# Fig S4) Regional IMP summary  #####
# x = fractional change
# y = regions
# split data points by crop
# panels dividied by phi

setwd(wddata)
fg2 = read.csv("summary.csv", header=TRUE,na.strings = c("#VALUE!", "#N/A", "N/A", "NA", ""))

# rename phi
fg2$phi =  as.factor(fg2$phi)
levels(fg2$phi) = c("Phi 2", "Phi 3", "Phi 4")

# change form long to wide (crteate columns for each phi)
fg2 = dcast(fg2, region + crop + phi ~ fact, value.var="value")

# fg2 = rowMeans(ALL[ ,c("IPM","MET","POP")],na.rm=TRUE) - not used

# leave IPM column
fg2=fg2[-which(names(fg2) %in% c("MET","POP"))]


# PLOT
library(ggplot2)
p = ggplot(fg2, aes(x = IPM, y = "region", color=crop)) 
p = p + geom_point() 
p = p + facet_wrap(~phi, ncol = 1,  scales="free")
p

# Save plot
# ppi = 300
# plotname = file.path(wdpng,paste("fig4",".png",sep = ""))
# png(filename=plotname,width=10*ppi, height=10*ppi, res=ppi )
# p
# dev.off()
# //////////////////









# ////////////////// Backyard: not worked through





### --------- FIGURE S2 : latitudinal median change in insect pest pressure at high (black), medium, (yellow) and low (blue) overwinter survival for maize (A), Rice (B), and Wheat (C).
# variables: IMP_xy


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
# Example of cicrle plot in ggplot:
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



# # functions
# source(file.path(wdfun,"color bar.R"))
# source(file.path(wdfun,"summarySE function.R"))
# source(file.path(wdfun,"multiplot.R"))




# Fig S5)  (Maize) projected future crop loss, current yield per hectare, and current yield loss due to pests, and projected future crop loss ----
# Variables
# Crop
# Yield loss
# Projected yield loss




# Fig S6)  (Rice) projected future crop loss, current yield per hectare, and current yield loss due to pests, and projected future crop loss ----
# Variables
# Crop
# Yield loss
# Projected yield loss




# Fig S7)  (Wheat) projected future crop loss, current yield per hectare, and current yield loss due to pests, and projected future crop loss ----
# Variables
# Crop
# Yield loss
# Projected yield loss




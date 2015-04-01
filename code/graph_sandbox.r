
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
# m1 = m1 + facet_wrap(~variable, ncol = 1,  scales="free")
m1 = m1 + facet_grid(variable~.)
m1 = m1 + guides(fill=guide_legend(title=NULL))

# add theme to map
m1 = m1 + theme_bw()
m1 = m1 + xlab("") + ylab("")

m1

# theme
m1 = m1 + theme_bw(base_family = "Helvetica")
m1 = m1 + theme(axis.title.x = element_text(face = "bold", color = "black", size = 12), axis.title.y = element_text(face = "bold", color = "black", size = 12), plot.title = element_text(face = "bold", color = "black", size = 16), strip.background = element_rect(fill = "white",colour = "white"), strip.text.x = element_text(colour = 'black', size = 12, face="bold"), strip.text.y = element_text(colour = 'black', size = 12, face="bold"))

m1 = m1 + theme(panel.grid = element_blank())

m1

# Save plot
# ppi = 300
# plotname = file.path(wdpng,paste("fig1",".png",sep = ""))
# png(filename=plotname,width=9*ppi, height=15*ppi, res=ppi )
# m1
# dev.off()
#####

#####
# --- Fig 1b) Summary of fractional change (horizontal plot, phi3 in the middle)

# data
fg = dat
fg = fg[,-1]

# sumarize
library(plyr)
fg = ddply(fg, c("region","crop","phi","fact"), summarise,value  = mean(value, na.rm=T)) # mean value / plot base

# re-arrange data: add columns for each phi
library(reshape2)
fg = dcast(fg, region + crop + fact ~ phi, value.var="value")
names(fg)[names(fg)=="2"] <- "phi2"
names(fg)[names(fg)=="3"] <- "phi3"
names(fg)[names(fg)=="4"] <- "phi4"

# melt phi for exploration
# pd = position_dodge(width = 1)
# ph = melt(fg, id.vars=c("region","crop","fact"))
# p = ggplot(ph,aes(y=value,x=crop,color=variable, group=variable)) + geom_point() + geom_point(position = pd) 
# p = p + ylim(ymin=-0.25, ymax=0.5)
# # p = p + coord_flip()
# p = p + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.x=element_line(colour="blue",size = 1))
# p + theme(panel.grid.major = element_blank())

# set facetting order
fg$fact = factor(fg$fact, levels=c('MET','POP','IPM'))

# change crop names
fg$crop = as.factor(fg$crop )
levels(fg$crop) = c("Maize", "Rice", "Wheat")

# //// PLOT
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

# init variables
pd = position_dodge(width = 1)

p = ggplot(fg, aes(x = region, y =phi3, color=crop)) # + geom_boxplot() 
p = p + geom_errorbar(aes(ymin=phi2, ymax=phi4), width=0.05, position=pd, width = 0.5) 
p = p + geom_point(position = pd)
p = p + coord_flip()
p = p + facet_grid(fact~.) 
p = p + ylim(xmin=-0.25, ymax=0.5)
p

# horizontal lines
p = p + geom_vline(xintercept=seq(0.5, length(unique(fg$region)), 1), lwd=0.2, colour="black")

# scales //  if colors exceed 11 use cols fun instead
n = 4
labs = levels(fg$crop)
p = p + scale_color_manual(name = "Crops", labels = labs, breaks=labs, values=brewer.pal("Spectral",n = n))

# theme
p = p + theme_bw(base_family = "Helvetica")
p = p + theme(axis.title.x = element_text(face = "bold", color = "black", size = 12), axis.title.y = element_text(face = "bold", color = "black", size = 12), plot.title = element_text(face = "bold", color = "black", size = 12), legend.position = c(1, 1), legend.justification = c(1, 1), strip.background = element_rect(fill = "white",colour = "white"), strip.text.x = element_text(colour = 'black', size = 12, face="bold"), strip.text.y = element_text(colour = 'black', size = 12, face="bold"))

p = p + theme(panel.grid = element_blank())

p
# Save plot
# ppi = 300
# plotname = file.path(wdpng,paste("fig1b_v1",".png",sep = ""))
# png(filename=plotname,width=8*ppi, height=12*ppi, res=ppi )
# p
# dev.off()

#####
#///// Combined plot p and m1

# remove title m1
m1 = m1 + theme(title = element_blank())

# remove xlab p
p = p + theme(axis.title.x = element_blank())

# place legends at the bottom (1a and 1b)
m1 = m1 + theme(legend.position="bottom")
p = p + theme(legend.position="bottom")

# remove facet names for 1b
p = p + theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())

# add fake tiles to 1b

m1

library(gridExtra)
grid.arrange(p,m1, ncol = 2, nrow = 1, widths = c(1,2), height = c(1, 1),main = textGrob("Fractional change in population metabolism",gp=gpar(fontsize=20,font=2)))









# --------------- 






p = ggplot(data=f1, aes(x = region, y = value, color=crop, group=phi)) 
# p = p + geom_line()
p = p + geom_line(stat='identity', position="dodge") 
p = p + coord_flip()
p = p + facet_grid(fact~.) 
# p = p + ylim(ymin=-0.5, ymax=0.5)
p



p = ggplot(data=f1, aes(x = region, y = value, color=crop, group=phi)) 
p = p +  geom_linerange(data = f2, mapping = aes(region, ymin = val_min, ymax = val_max, color = crop))
p = p + coord_flip()
p = p + facet_grid(fact~.) 
p = p + ylim(ymin=-0.5, ymax=0.5)
p




+ geom_line() 

aes(x = Jahr, y = Wert, color = Altersgr, group = Altersgr)) + 
  geom_point() + geom_line() + 
  
  
  + geom_line(position = pd)

p = p + geom_linerange(data = f2, aes(region, ymin = val_min, ymax = val_max, color = crop),position = pd)
p

p + geom_linerange(aes(as.factor(region),ymin=min(value),ymax=value,color=crop))

+ geom_bar(position = pd) 
p = p + coord_flip()
p = p + facet_grid(fact~.) 
p = p + ylim(ymin=-0.5, ymax=0.5)
p = p + xlab(label = "") + ylab(label = "Fractional change") 
p



-----
  # melt phi for exploration
  # pd = position_dodge(width = 1)
  # ph = melt(fg1b, id.vars=c("region","crop","fact"))
  # p = ggplot(ph,aes(y=value,x=crop,color=variable, group=variable)) + geom_point() + geom_point(position = pd) 
  # p = p + ylim(ymin=-0.25, ymax=0.5)
  # # p = p + coord_flip()
  # p = p + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.x=element_line(colour="blue",size = 1))
  # p + theme(panel.grid.major = element_blank())
  ---

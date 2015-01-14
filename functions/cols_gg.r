# This code chunk sets the color palette and custom gg variables

# Create custom ramp palette based on "spectral"
library(RColorBrewer)
pal = brewer.pal("Spectral",n = 11)
cols = colorRampPalette(pal)

# custom palette for crops
palcrop = c("#9E0142","#F88D51","#5E4FA2")

# HOW TO USE:
# editing gg scales [if colors exceed 11 use cols() fun instead]
# labs = levels(fg$crop)
# p = p + scale_color_manual(name = "Crops", labels = labs, breaks=labs, values=brewer.pal("Spectral",n = n))

# theme
library(ggplot2)
mygg = theme_bw(base_family = "Helvetica") + theme(axis.title.x = element_text(face = "bold", color = "black", size = 12), axis.title.y = element_text(face = "bold", color = "black", size = 12), plot.title = element_text(face = "bold", color = "black", size = 14), strip.background = element_rect(fill = "white",colour = "white"), strip.text.x = element_text(colour = 'black', size = 12, face="bold"), strip.text.y = element_text(colour = 'black', size = 12, face="bold"))

# p = p + theme(panel.grid = element_blank()) # edit grids

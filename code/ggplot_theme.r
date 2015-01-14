# Prepare color palette
library(RColorBrewer)
pal <- brewer.pal("Spectral",n = 11)
cols<- colorRampPalette(pal)

# scales //  if colors exceed 11 use cols fun instead
n = 4
labs = levels(fg$crop)
p = p + scale_color_manual(name = "Crops", labels = labs, breaks=labs, values=brewer.pal("Spectral",n = n))

# theme
p = p + theme_bw(base_family = "Helvetica")
p = p + theme(axis.title.x = element_text(face = "bold", color = "black", size = 12), axis.title.y = element_text(face = "bold", color = "black", size = 12), plot.title = element_text(face = "bold", color = "black", size = 12), legend.position = c(1, 1), legend.justification = c(1, 1), strip.background = element_rect(fill = "white",colour = "white"), strip.text.x = element_text(colour = 'black', size = 12, face="bold"), strip.text.y = element_text(colour = 'black', size = 12, face="bold"))

# p = p + theme(panel.grid = element_blank()) # edit grids

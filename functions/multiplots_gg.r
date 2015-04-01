library (ggplot2)
# placeholder plot - prints nothing at all (if needed)
empty_space <- ggplot() + geom_point(aes(1, 1), colour = "white") + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())

# library(gridExtra)

# EXAMPLE
# grid.arrange(plot_top, empty, scatter, plot_right, ncol = 2, nrow = 2, widths = c(4, 1), height = c(1, 4))

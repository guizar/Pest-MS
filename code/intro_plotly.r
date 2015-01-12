# Introducing plotly

# First update packages
update.packages()

# install devtools  (see if it's instaled already)
install.packages("devtools")  # so we can install from github

library("devtools")
install_github("ropensci/plotly")  # plotly is part of ropensci

library(plotly)

py <- plotly(username="guizar", key="v0r8an3h4p")  # open plotly connection

py$ggplotly()

m1
# R > 4.4
install.packages("Matrix")
install.packages("igraph")
install.packages("remotes")
install.packages("foreach")
install.packages("doParallel")
install.packages("dplyr")
remotes::install_github("https://github.com/asgr/celestial")
remotes::install_github("https://github.com/asgr/cmaeshpc")
remotes::install_github("https://github.com/asgr/Highlander")
Sys.setenv(NOT_CRAN = "true")
install.packages("arrow")
remotes::install_local("./FoFR")

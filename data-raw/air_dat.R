## code to prepare `air_dat` dataset goes here
library(CARBayesdata)
library(spdep)
data("GGHB.IG", package = "CARBayesdata")
data("pollutionhealthdata", package = "CARBayesdata")
raw_airdat <- list(spatial = GGHB.IG, df = pollutionhealthdata)

usethis::use_data(raw_airdat)

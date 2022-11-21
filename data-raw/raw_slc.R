## code to prepare `raw_slc` dataset goes here
library(CARBayesdata)
library(shapefiles)
library(sp)
data(lipdata)
data(lipdbf)
data(lipshp)

library(CARBayes)
lipdb <- lipdbf
lipdb$dbf <- lipdbf$dbf[ ,c(2,1)]
data.combined <- combine.data.shapefile(data=lipdata, shp=lipshp, dbf=lipdb)

W.nb <- poly2nb(data.combined, row.names = rownames(lipdata))
W.mat <- nb2mat(W.nb, style="B")
library(spdep)
W.nb <- poly2nb(propertydata.spatial, row.names = rownames(df))
W <- nb2mat(W.nb, style="B")
raw_slcdat <- list(df = df, W = W, y = df$logprice,
                   formula = logprice~ns(crime, 3) + rooms + sales + factor(type) + logdriveshop)
usethis::use_data(raw_slcdat, overwrite = TRUE)

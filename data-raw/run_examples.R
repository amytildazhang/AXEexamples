# Code to run AXE and other LCO methods for the examples in the paper are here.
# Comments are provided for re-creating the eight schools example;
# other examples follow the same structure.


for (file in list.files("R")) {
  source(file.path("R", file))
}


###
#  Eight schools example
###

# accesses and/or saves data for the example + scenario parameters
# e.g. data scaling factor alpha or number of clusters
eight <- prep_eight()

# runs MCV and saves estimates; also records the time
# examine using str(eight$cv_yhats)
eight$cv_yhats <- mcv_eight()

# fits models to full data and saves posterior means of variance parameters
eight$posteriors <- pfit_eight()

# runs AXE and saves estimates; also records the time
eight$axe_yhats <- axe_eight()
# usethis::use_data(eight, overwrite = T) # saves the object `eight` for use within the package


###
#  Radon example (
###
radon_1 <- prep_radon_full()
radon_1$cv_yhats <- mcv_radon_full()
radon_1$posteriors <- pfit_radon_full()
radon_1$axe_yhats <- axe_radon_full()
# usethis::use_data(radon_1, overwrite = T)


###
#  Radon subsets
###

radon_2 <- prep_radon_simul()

# this is structured slightly differently, just due to the sheer amount of time it takes to run
# First, get posterior samples for all subsets and save in individual files (not included, requires 1.18Gb)
radon_2$posterior_summary <- pfit_radon_simul(use_saved = F)
# Next, run each LCO method individually on posterior samples
radon_2$lco <- list()
radon_2$lco$ghst_yhats <- lco_radon_simul(method = "ghost")
radon_2$lco$iisfm_yhats <- lco_radon_simul(method = "iis_fm")
radon_2$lco$iisim_yhats <- lco_radon_simul(method = "iis_im")

# re-combine LCO methods with posteriors, to be in line with other datasets
radon_2$posteriors <- dplyr::full_join(
  radon_2$posterior_summary,
  radon_2$lco$iisfm_yhats
) %>%
  dplyr::full_join(radon_2$lco$iisim_yhats) %>%
  dplyr::full_join(radon_2$lco$ghst_yhats)

radon_2$cv_yhats <- mcv_radon_simul()

# no real change found from comparing to IIS variance or MAP
# # need radon_2$posteriors before running this
# radon_2$axe_yhats <- dplyr::full_join(
#     axe_radon_simul(var = "iis") %>%
#         dplyr::select(perc, iter, n_clusters, model, idx, yhat_axe) %>%
#         dplyr::rename(yhat_axe_iis = yhat_axe),
#     axe_radon_simul(var = "map") %>%
#         dplyr::select(perc, iter, n_clusters, model, idx, yhat_axe) %>%
#         dplyr::rename(yhat_axe_map = yhat_axe)
# ) %>% dplyr::full_join(axe_radon_simul(var = "axe"))

radon_2$axe_yhats <- (axe_radon_simul(var = "axe"))
# usethis::use_data(radon_2, overwrite = T)




# ESP
lol <- prep_lol()
lol$posteriors <- pfit_lol()
lol$cv_yhats <- mcv_lol()
lol$axe_yhats <- axe_lol()
# usethis::use_data(lol, overwrite = T)


# SLC
slc <- list()
slc$data <- prep_slc()
slc$cv_yhats <- mcv_slc(n_cores = 1)
slc$posteriors <- pfit_slc()
slc$axe_yhats <- axe_slc()
# usethis::use_data(slc, overwrite = T)


# SRD
air <- list()
air$data <- prep_air()
air$cv_yhats <- mcv_air()
air$posteriors <- pfit_air()
air$axe_yhats <- axe_air()
# usethis::use_data(air, overwrite = T)

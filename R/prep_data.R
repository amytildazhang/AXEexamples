# functions to prep data for eight schools example


# Supporting functions for prepping data


create_schema <- function(c_county, c_size, dperc, nclusters, sizes,
                          seed, max_combns, max_times = 1, ...) {
    set.seed(seed)
    purrr::map(dperc, function(perc) {
        print(perc)
        total <- round(c_size/perc)


        purrr::map(nclusters, function(nc) {
            print(nc)

            n <- nc - 1
            combns <- RcppAlgos::comboGeneral(sizes, m = n, upper = 1e3,
                                              limitConstraints = c(round(0.9 * total), round(1.1 * total)),
                                              constraintFun = "sum",
                                              comparisonFun = c(">=", "<="))

            if (nrow(combns) == 0) {
                return(NULL)
            }


            if (nrow(combns) > max_combns) {
                # sample at most 60 rows
                vals <- rowSums(combns)
                abs_dif <- abs(vals - total)
                prob_weights <- 1/abs_dif
                prob_weights[is.infinite(prob_weights)] <- 1


                chosen <- sample(1:nrow(combns), max_combns, prob = prob_weights)

                combns <- combns[chosen,, drop = F]
            }


            # for each combination, choose counties
            do.call(
                "rbind",
                purrr::map(1:nrow(combns), function(i) {
                    nums_chosen <- combns[i, ]
                    names_chosen <-  purrr::map(unique(nums_chosen), function(size) {
                        n_size <- sum(nums_chosen == size)
                        options <- names(sizes)[sizes == size]
                        if (length(options) >= n_size) {
                            t(combn(options, n_size))
                        } else {
                            NA
                        }

                    })


                    lengths <- purrr::map_int(names_chosen, nrow)
                    total_length <- prod(lengths)


                    county_combn <- do.call("cbind", purrr::map(1:length(names_chosen), function(j) {
                        if (j > 1) {
                            n_each <- prod(lengths[1:(j - 1)])
                        } else {
                            n_each <- 1
                        }

                        n_rep <- total_length/(n_each * lengths[j])
                        if (round(n_rep) != n_rep) stop("Unequal divisions")
                        idx <- rep(1:nrow(names_chosen[[j]]), times = n_rep, each = n_each)
                        names_chosen[[j]][idx, , drop = F]
                        # rep_mat(names_chosen[[j]], times = n_rep, each = n_each)

                    }))

                    if (ncol(county_combn) != n) stop(sprintf("Incorrect number of columns in %s", i))

                    has_na <-  apply(county_combn, 1, function(row) all(is.na(row)))
                    county_combn <- county_combn[!has_na, , drop = F]


                    # sample 1
                    county_combn[sample(1:nrow(county_combn), max_times), , drop = F]
                })
            )

        }) %>% purrr::set_names(sprintf("cluster%s", nclusters))
    }) %>% purrr::set_names(sprintf("perc%s", dperc))

}

####### final functions to prep


#' Data used in paper
#'
#' `prep_*()` returns a list with objects used to fit posteriors and obtain LCO
#'      approximations. See README.
#'
#' @param data_scale Eight schools. Value of alpha (see paper)
#' @param c_county Radon subsets. Name of test county
#' @param source Radon subsets. List with dataframe to use (default is `radon_1`).
#' @param data_perc Radon subsets. Test data percentage.
#' @param n_clusters Radon subsets. Number of counties in data subset.
#' @param max_combn Radon subsets. Max number of subsets for each `data_perc` and `n_clusters.`
#' @param max_times Radon subsets. Number of times to iterate.
#' @param seed Radon subsets. Seed for subset generation.
prep_eight <- function(data_scale = seq(0.1, 4, by = 0.1)) {
    J <- 8
    school_idx <- 1:J

    dat <- data.frame(
        school = factor(letters[school_idx]),
        y = c(28,  8, -3,  7, -1,  1, 18, 12),
        sd = c(15, 10, 16, 11,  9, 11, 10, 18)
    )
    X <- model.matrix(~school, data = dat,
                      contrasts.arg = list(school = contrasts(dat$school, contrasts = F)))
    schools_dat <- list(
        J = J,
        y = dat$y,
        school_idx = school_idx,
        sigma = dat$sd,
        sigma_idx = 1:J,
        N = nrow(X),
        H = 0,
        train_idx = 1:nrow(dat),
        test_idx = vector()
    )



    list(df = dat, stan_list = schools_dat, X = X, data_scale = data_scale,
         school_idx = school_idx)
}



#' @describeIn prep_eight Function to get data for Radon example.
prep_radon_full <- function() {
    raw_radon <- list()
    raw_radon$data <- rstanarm::radon

    raw_radon$models_mcv <- list(
        log_radon ~ 1 + (1 | county),
        log_radon ~ floor + (1 | county),
        log_radon ~ floor + log_uranium + (1 | county)
    )

    raw_radon$models_axe <- list(
        ~ county,
        ~ -1 + floor +  county,
        ~ -1 + floor + log_uranium +  county
    )



    raw_radon
}



#' @describeIn prep_eight Function to get data for Radon subsets example.
prep_radon_simul <- function(c_county = "OLMSTED", source = radon_1,
                             data_perc =  c(0.3, 0.4, 0.5, 0.6, 0.7),
                             n_clusters = c(3, 4, 6, 9, 12),
                             max_combn = 60, max_times = 1, seed = 9871) {
    clusters <- source$data$county

    c_size <- sum(source$data$county == c_county)
    totals <- c_size/data_perc
    c_size/totals

    sizes <- table(as.character(clusters[clusters != c_county]))

    schema <- create_schema(c_county, c_size, data_perc, n_clusters,
                            sizes, seed, max_combn, max_times)

    list(
        schema = schema, c_county = c_county, c_size = c_size,
        data_perc = data_perc, n_clusters = n_clusters,
        data = source$data, models_mcv = source$models_mcv,
        models_axe = source$models_axe
    )
}



#' @describeIn prep_eight Function to get data for Esports example.
prep_lol <- function() {
    df <- read.csv("data-raw/lolesp.csv") %>%
        dplyr::filter(league == "LCS", position != "team") %>%
        dplyr::select(kills, team, player, champion, gamelength, dpm, wpm, earned.gpm, teamkills, position) %>%
        dplyr::mutate_at(dplyr::vars(player, champion, position, team), function(vec) {
            factor(stringr::str_replace_all(vec, "\\W", ""))
        }) %>%
        dplyr::mutate(teamkills = teamkills - kills,
                      log_dpm = log(dpm),
                      log_egpm = log(earned.gpm),
                      idx = 1:dplyr::n())


    contr <- contrasts_for_pooling(df, c("player", "champion"))
    X <- model.matrix(dpm ~  position  + team +
                          log_dpm + log_egpm  + champion + player, data = df, contrasts = contr)


    list(data = df, X = X, loops = as.character(df$player))

}


#
# combns <- function(names) {
#     combn(sort(names), 2)
# }


#' @describeIn prep_eight Function to get data for Scottish Lung Cancer (SLC) example.
prep_slc <- function() {
    c(raw_slcdat, list(loops = 1:data$n))
}



#' @describeIn prep_eight Function to get data for Scottish Respiratory Disease example.
prep_air <- function() {
    df <- raw_airdat$df %>% dplyr::mutate(
        SMR = observed / expected,
        logSMR = log(SMR)
    )

    SMR.av <- dplyr::group_by(df, IG) %>%
        dplyr::summarise(SMR.mean = mean(SMR), .groups = "drop")
    raw_airdat$spatial@data$SMR <- SMR.av$SMR.mean

    W.nb <- spdep::poly2nb(air_dat$spatial, row.names = SMR.av$IG)
    W.list <- spdep::nb2listw(W.nb, style = "B")
    W <- spdep::nb2mat(W.nb, style = "B")
    formula <- observed ~ offset(log(expected)) + jsa + price + pm10

    X <- model.matrix(~jsa + price + pm10, data = df)
    list(
        W = W, log_offset = log(df$expected), df = df,
        loops = df$IG,  X = X, formula = formula
    )

}

# These are functions which produce graphs and tables used in the paper

#' Rename character vector with methods
#'
#' Change internal method names to plot output
#'
#' @export
rename_methods <- function(vec) {
  vec <- stringr::str_replace(vec, "_fm$", "") |>
    stringr::str_replace("ghst_c", "GHOST") |>
    stringr::str_to_upper() |>
    stringr::str_replace("_", "-") |>
    stringr::str_replace("^IIS", "iIS") |>
    stringr::str_replace("VT", "VEHTARI")


  dplyr::case_when(
    stringr::str_detect(vec, "PSIIS") ~
      stringr::str_replace(vec, "PSIIS", "iIS"),
    vec == "PSIS-C" ~ "PSIS-LOO",
    vec == "CV-SIJ" ~ "PPMode",
    TRUE ~ vec
  )
}

make_results_df <- function(posteriors, axe, mcv_vals) {
  res <- dplyr::full_join(
    posteriors$df,
    axe$df,
    by = c()
  ) |>
    dplyr::full_join(mcv_vals)

  # add AXE time
  J <- length(unique(res$loop))
  res |>
    dplyr::group_by(loop) |>
    dplyr::mutate(time_axe = axe$time / J)
}

# time for MCV
time_eight <- function(info = eight) {
  times_df <- purrr::map_df(seq_along(info$cv_yhats), function(i) {
    data.frame(
      idx = info$cv_yhats[[i]]$idx,
      time_cv = info$cv_yhats[[i]]$time_cv,
      data_scale = info$cv_yhats[[i]]$data_scale,
      loop = 1:8,
      time_iis_fm = info$posteriors[[i]]$time_iis_fm,
      time_iis_im = info$posteriors[[i]]$time_iis_im,
      time_ghst_c = info$posteriors[[i]]$time_iis_im
    )
  }) |>
    dplyr::select(loop, dplyr::starts_with("time")) |>
    unique() |>
    dplyr::summarise_all(sum) |>
    dplyr::select(-loop)


  times_df |> dplyr::mutate(time_axe = info$axe_yhats$time)
}

#' Dataframe collating LCO, MCV, and AXE results
#'
#' `results_*()` combines the output from `mcv_*`, `axe_*`, and `pfit_*`
#'     into one dataframe.
#' @return data.frame
#' @export
results_eight <- function(info = eight, mcv_vals = eight$cv_yhats,
                          axe_vals = eight$axe_yhats,
                          posteriors = eight$posteriors) {

  dplyr::full_join(
    posteriors |> dplyr::rename(idx = i),
    axe_vals |> dplyr::rename(yhat_axe = axe_yhat)
  ) |>
    dplyr::full_join(
      dplyr::select(mcv_vals, cv_yhat, time_cv,
                    yhat_cv_sij, idx, data_scale) |>
        dplyr::rename(yhat_cv = cv_yhat)
    ) |>
    dplyr::mutate(
      school_idx = letters[school_idx],
      y = info$df$y[loop] * data_scale
    )
}


#' @describeIn results_eight
#' @export
results_radon_full <- function(info = radon_1,
                               posteriors = radon_1$posteriors,
                               mcv_vals = radon_1$cv_yhats,
                               axe_vals = radon_1$axe_yhats) {
  purrr::map_df(seq_along(info$models_axe), function(i) {
    resdf <- dplyr::full_join(
      info$data |> dplyr::mutate(idx = 1:nrow(info$data)),
      mcv_vals[[i]] |> dplyr::select(-loop),
      by = c()
    ) |>
      dplyr::full_join(
        posteriors[[i]] |> dplyr::mutate(loop = as.numeric(loop)),
        by = c("idx" = "i")
      )
    J <- length(unique(resdf$loop))

    resdf |> dplyr::mutate(time_axe = axe_vals$time / J)
  }) |>
    dplyr::full_join(axe_vals$axe, by = c("idx" = "i", "model")) |>
    dplyr::rename(y = log_radon)
}

time_radon_full <- function(info = radon_1) {
  cvdf <- dplyr::bind_rows(info$cv_yhats)
  posdf <- dplyr::bind_rows(info$posteriors) |>
    dplyr::mutate(loop = as.numeric(loop))

  dplyr::full_join(cvdf, posdf) |>
    dplyr::group_by(loop, model) |>
    dplyr::select(dplyr::starts_with("time")) |>
    unique() |>
    dplyr::ungroup() |>
    dplyr::summarise_at(dplyr::vars(dplyr::starts_with("time")), sum) |>
    dplyr::mutate(time_axe = info$axe_yhats$time)
}

#' @describeIn results_eight
#' @export
results_radon_simul <- function(info = radon_2, mcv_vals = radon_2$cv_yhats,
                                axe_vals = radon_2$axe_yhats,
                                posteriors = radon_2$posteriors) {
  agg <- dplyr::full_join(
    mcv_vals,
    axe_vals,
    by = c()
  ) |>
    dplyr::full_join(posteriors) |>
    dplyr::rename(elpd_cv_rst = rst_elpd) |>
    dplyr::group_by(iter, perc, n_clusters, model, n_tot)


  agg
}

time_radon_simul <- function(info = radon_2) {
  dplyr::full_join(info$cv_yhats, info$posteriors) |>
    dplyr::full_join(info$axe_yhats) |>
    dplyr::group_by(iter, perc, n_clusters, model) |>
    dplyr::select(dplyr::contains("time")) |>
    unique() |>
    dplyr::ungroup() |>
    dplyr::summarise_at(dplyr::vars(dplyr::contains("time")), sum)
}


#' @describeIn results_eight
#' @export
results_lol <- function(info = lol, posteriors = lol$posteriors,
                        axe = lol$axe_yhats, mcv_vals = lol$cv_yhats) {
  mcv <- do.call("rbind", mcv_vals)
  make_results_df(posteriors, axe, mcv) |>
    dplyr::full_join(info$data |>
                       dplyr::mutate(idx = 1:dplyr::n()) |>
                       dplyr::select(idx, kills) |>
                       dplyr::rename(y = kills))
}


#' @describeIn results_eight
#' @export
results_slc <- function(info = slc, posteriors = slc$posteriors,
                        axe = slc$axe_yhats, mcv_vals = slc$cv_yhats) {
  make_results_df(posteriors, axe, mcv_vals) |>
    dplyr::full_join(
      data.frame(
        y = info$data$y,
        idx = seq_along(info$data$y)
      )
    )
}



time_slc <- function(info = slc) {
  dplyr::full_join(info$cv_yhats, info$posteriors$df) |>
    dplyr::group_by(loop) |>
    dplyr::select(dplyr::contains("time")) |>
    unique() |>
    dplyr::ungroup() |>
    dplyr::select(-loop) |>
    dplyr::summarise_all(sum) |>
    dplyr::mutate(time_axe = info$axe_yhats$time)
}



#' @describeIn results_eight
#' @export
results_air <- function(info = air, posteriors = air$posteriors,
                        axe = air$axe_yhats, mcv_vals = air$cv_yhats) {
  make_results_df(posteriors, axe, mcv_vals) |>
    dplyr::full_join(
      data.frame(
        y = info$data$df$observed,
        idx = 1:nrow(info$data$df)
      ),
      by = "idx"
    )
}


compare_time <- function(resdf) {
  resdf |>
    dplyr::group_by(loop) |>
    dplyr::select(dplyr::contains("time")) |>
    unique() |>
    dplyr::summarise_at(-loop, sum)
}

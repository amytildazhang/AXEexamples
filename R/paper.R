
# These are functions which produce graphs and tables used int he paper

#' Rename character vector with methods
#'
#' Change internal method names to plot output
#'
#' @export
rename_methods <- function(vec) {
  vec <- stringr::str_replace(vec, "_fm$", "-C") %>%
    stringr::str_replace("_im$", "-A") %>%
    stringr::str_replace("ghst_c", "GHOST") %>%
    stringr::str_replace("_", "-") %>%
    stringr::str_to_upper()


  dplyr::case_when(
    stringr::str_detect(vec, "WAIC") | stringr::str_detect(vec, "DIC") ~ sprintf("i%s", vec),
    stringr::str_detect(vec, "PSIIS") ~ stringr::str_replace(vec, "PSIIS", "iIS"),
    vec == "PSIS-C" ~ "PSIS-LOO",
    T ~ vec
  )
}

#' Colors for each method in plots
#'
#' Returns named vector of colors for each method, for use in `ggplot2::scale_color_manual()`.`
#'
#' @export
method_colors <- function(type = c("elpds", "yhats"),
                          methods = c(
                            "AXE", "AXE-IIS", "AXE-MAP", "GHOST",
                            "iIS-C", "iIS-A",
                            "DIC-C", "DIC-A", "WAIC-C", "WAIC-A"
                          )) {
  if (type == "elpd") {
    cols <- RColorBrewer::brewer.pal(length(methods), "Paired")
  } else if (type == "yhat") {
    cols <- RColorBrewer::brewer.pal(length(methods), "Dark2")
  }


  names(cols) <- unique(methods)
  cols
}


#' Create a dataframe which compares methods.
#'
#' Creates dataframe with LRRs, as defined in paper.
#'
#' @export
df_compare_methods <- function(df, type = c("yhat", "elpd")) {
  # df requirements: must have loop column, actual response value y, and
  # columns with method approximations must contain "yhat" or "elpd
  if (type == "yhat") {
    df %>%
      # calculate RMSE for each loop
      dplyr::summarise_at(
        dplyr::vars(dplyr::contains("yhat")),
        ~ sqrt(mean((. - y)^2))
      ) %>%
      # set aside ground truth MCV values
      dplyr::rename(ref_y = yhat_cv) %>%
      tidyr::pivot_longer(
        cols = dplyr::contains("yhat"),
        names_to = "method", values_to = "yhat"
      ) %>%
      dplyr::mutate(
        dif = (log((yhat / ref_y))),
        method = stringr::str_replace(method, "yhat_", "") %>%
          rename_methods()
      )
  } else {
    df %>%
      dplyr::summarise_at(
        dplyr::vars(dplyr::contains("elpd")),
        function(elpd) {
          if (length(unique(elpd)) == 1) {
            elpd[1]
          } else {
            sum(elpd)
          }
        }
      ) %>%
      dplyr::rename(ref_elpd = elpd_cv) %>%
      tidyr::pivot_longer(
        cols = dplyr::starts_with("elpd"),
        names_to = "method", values_to = "elpd"
      ) %>%
      dplyr::mutate(
        dif = abs((ref_elpd - elpd) / ref_elpd),
        method = stringr::str_replace(method, "elpd_", "") %>%
          rename_methods()
      )
  }
}

#' Plot comparison of methods (LRR)
#'
#' Assumes use of `df_compare_methods`
#'
#' @param df Output from `df_compare_methods`
#' @export
plot_compare_methods <- function(df, type = c("yhat", "elpd")) {
  if (type == "yhat") {
    colors <- method_colors(
      type = "yhat",
      methods = c("AXE", "GHOST", "iIS-C", "PSiIS-C", "iIS-A", "PSiIS-A")
    )
    ylab <- "log | Y-hat Method/MCV|"
  } else {
    colors <- method_colors(
      type = "elpd",
      methods = c(
        "iDIC-C", "iDIC-A",
        "iWAIC-C", "iWAIC-A",
        "iIS-C", "iIS-A",
        "GHOST"
      )
    )
    ylab <- "|(ELPD Method - MCV)/MCV|"
  }

  ggplot2::ggplot(df, ggplot2::aes(x = method, y = dif)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      # text = ggplot2::element_text(family = "CMU Serif"),
      legend.key.size = ggplot2::unit(0.3, "lines"),
      legend.background = ggplot2::element_rect(
        fill = ggplot2::alpha("white", 0)
      )
    ) +
    ggplot2::labs(
      y = ylab,
      x = NULL,
      fill = NULL, color = NULL
    )
}

time_eight <- function(info = eight) {
  times_df <- purrr::map_df(1:length(info$cv_yhats), function(i) {
    data.frame(
      idx = info$cv_yhats[[i]]$idx,
      time_cv = info$cv_yhats[[i]]$time_cv,
      data_scale = info$cv_yhats[[i]]$data_scale,
      loop = 1:8,
      time_iis_fm = info$posteriors[[i]]$time_iis_fm,
      time_iis_im = info$posteriors[[i]]$time_iis_im,
      time_ghst_c = info$posteriors[[i]]$time_iis_im
    )
  }) %>%
    dplyr::select(loop, dplyr::starts_with("time")) %>%
    unique() %>%
    dplyr::summarise_all(sum) %>%
    dplyr::select(-loop)


  times_df %>% dplyr::mutate(time_axe = info$axe_yhats$time)
}

#' Dataframe collating LCO, MCV, and AXE results
#'
#' `results_*()` combines the output from `mcv_*`, `axe_*`, and `pfit_*`
#'     into one dataframe.
#' @return data.frame
#' @export
results_eight <- function(info = eight, mcv_vals = eight$cv_yhats, axe_vals = eight$axe_yhats$axe,
                          posteriors = eight$posteriors) {
  purrr::map_df(1:length(axe_vals), function(i) {

    # sanity checks
    stopifnot(axe_vals[[i]]$data_scale == mcv_vals[[i]]$data_scale)
    stopifnot(posteriors[[i]]$data_scale[1] == mcv_vals[[i]]$data_scale)


    # convenience objects
    cvhat <- mcv_vals[[i]]
    axehat <- axe_vals[[i]]
    psishat <- posteriors[[i]]

    # save estimates for Xbeta via all three methods
    posteriors[[i]] %>%
      dplyr::mutate(
        school_idx = letters[school_idx],

        yhat_cv = cvhat$cv_yhat[cvhat$idx],
        yhat_axe = axehat$axe_yhat[axehat$idx],
        elpd_mcv_c = cvhat$cv_elppd,
        elpd_mcv_m = cvhat$cv_lco,
        time_cv = cvhat$time_cv[cvhat$idx]
      )
  }) %>%
    dplyr::mutate(y = info$df$y[loop] * data_scale)
}

#
# figure_eight <- function(results) {
#
#
#     p <-
#         results %>%
#         # reshape the table and re-label categories for plotting
#         tidyr::gather(key = "yhat", value = "val", dplyr::starts_with("yhat_")) %>%
#         dplyr::mutate(
#             yhat = dplyr::case_when(yhat == "yhat_axe" ~ "AXE",
#                                     yhat == "yhat_cv" ~ "MCV",
#                                     yhat == "yhat_psis" ~ "PSIS_c",
#                                     yhat == "yhat_psis_m" ~ "PSIS_m",
#                                     yhat == "yhat_is" ~ "IS",
#                                     yhat == "yhat_ghst" ~ "GHOST") %>%
#                 factor(levels = c("MCV", "PSIS_c", "PSIS_m", "IS", "GHOST", "AXE"), ordered = T)
#         ) %>%
#         # plot base
#         ggplot2::ggplot(
#             ggplot2::aes(
#                 x = data_scale, y = abs(val),
#                 color = yhat, shape = yhat)
#         ) +
#         ggplot2::geom_point(alpha = 0.7) +
#         ggplot2::geom_line() +
#
#         # colors and shapes
#         ggplot2::scale_color_manual(values = c("black", RColorBrewer::brewer.pal(6, "Dark2") )) +
#         ggplot2::scale_shape_manual(values = c(16, rep(1, 6))) +
#
#         lemon::facet_rep_wrap(~letters[1:8], ncol = 5) +
#
#         # labels and theming
#         ggplot2::labs(
#             color = NULL, shape = NULL,
#             x = expression(alpha),
#             y = expression(hat(Y)[pp])
#         ) +
#         ggplot2::theme_bw() +
#         ggplot2::theme(
#             panel.border = ggplot2::element_blank(),
#             axis.line = ggplot2::element_line()
#         )
#
#
#     # move legend into empty panel
#     lemon::reposition_legend(p, position = "bottom right", panel = "panel-3-2",
#                              offset = 0.15)
#
#
# }
#


#########

#' @describeIn results_eight
#' @export
results_radon_full <- function(info = radon_1,
                               posteriors = radon_1$posteriors,
                               mcv_vals = radon_1$cv_yhats,
                               axe_vals = radon_1$axe_yhats) {
  purrr::map_df(1:length(info$models_axe), function(i) {
    resdf <- dplyr::full_join(
      info$data %>% dplyr::mutate(idx = 1:nrow(info$data)),
      mcv_vals[[i]] %>% dplyr::select(-loop),
      by = c()
    ) %>%
      dplyr::full_join(
        posteriors[[i]] %>% dplyr::mutate(loop = as.numeric(loop)),
        by = c("idx" = "i")
      )
    J <- length(unique(resdf$loop))

    resdf %>% dplyr::mutate(time_axe = axe_vals$time / J)
  }) %>%
    dplyr::full_join(axe_vals$axe, by = c("idx" = "i", "model")) %>%
    dplyr::rename(y = log_radon)
}

time_radon_full <- function(info = radon_1) {
  cvdf <- dplyr::bind_rows(info$cv_yhats)
  posdf <- dplyr::bind_rows(info$posteriors) %>% dplyr::mutate(loop = as.numeric(loop))

  dplyr::full_join(cvdf, posdf) %>%
    dplyr::group_by(loop, model) %>%
    dplyr::select(dplyr::starts_with("time")) %>%
    unique() %>%
    dplyr::ungroup() %>%
    dplyr::summarise_at(dplyr::vars(dplyr::starts_with("time")), sum) %>%
    dplyr::mutate(time_axe = info$axe_yhats$time)
}

# figure_radon_full <- function(results) {
#  p_yhat <-  results %>%
#         dplyr::select(-yhat_dif) %>%
#         dplyr::rename(ref_val = yhat_cv) %>%
#         tidyr::gather("yhat", "value", tidyselect::starts_with("yhat")) %>%
#         dplyr::group_by(county, model, yhat) %>%
#         dplyr::summarise(rmse = sqrt(mean((value - ref_val)^2))) %>%
#         dplyr::mutate(
#             yhat = stringr::str_replace(yhat, "yhat_", "") %>%
#                 rename_methods(),
#             model = sprintf("Model %s", model)
#         ) %>%
#
#         ggplot2::ggplot() +
#         ggplot2::geom_point(ggplot2::aes(x = county, y = rmse, color = yhat)) +
#         ggplot2::facet_wrap(~model) +
#         ggplot2::theme_bw() +
#         ggplot2::theme(axis.text.x = ggplot2::element_blank(),
#               axis.ticks = ggplot2::element_blank()) +
#         ggplot2::geom_hline(yintercept = 0, color = 'darkgray') +
#        ggplot2::labs(color = NULL, y = "RMSE from MCV, by county", x = "County")
#
#    # draw with grid::grid.draw
#    p_yhat <- lemon::reposition_legend(p_yhat, x = 0.95, y = 0.95, just = c(1, 1), panel = "panel-1-1")
#
#
#
#
#    p_elpd <-
#        results %>%
#            dplyr::select(-elpd_iis_m, -elpd_iis_c) %>%
#            dplyr::rename(ref_val = elpd_cv) %>%
#        tidyr::gather("elpd", "value", tidyselect::starts_with("elpd_")) %>%
#        dplyr::group_by(county, model, elpd) %>%
#        dplyr::summarise(
#            cv_elpd = sum(ref_val),
#            elpd_dif = sum(value) - sum(ref_val)
#        ) %>%
#        dplyr::mutate(
#            elpd = stringr::str_replace(elpd, "elpd_", "") %>%
#                rename_methods(),
#            model = sprintf("Model %s", model)
#        ) %>%
#        ggplot2::ggplot() +
#        ggplot2::geom_point(ggplot2::aes(x = county, y = elpd_dif, color = elpd)) +
#        ggplot2::facet_wrap(~model) +
#        ggplot2::theme_bw() +
#        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
#                       axis.ticks = ggplot2::element_blank()) +
#        ggplot2::geom_hline(yintercept = 0, color = 'darkgray') +
#        ggplot2::labs(color = NULL, y = "Difference from MCV ELPD, by county", x = "County")
#    p_elpd <- lemon::reposition_legend(p_elpd, x = 0.95, y = 0.95, just = c(1, 1), panel = "panel-1-1")
#
#
#    list(p_yhat = p_yhat, p_elpd = p_elpd)
#
# }
#

#####################

#' @describeIn results_eight
#' @export
results_radon_simul <- function(info = radon_2, mcv_vals = radon_2$cv_yhats,
                                axe_vals = radon_2$axe_yhats, posteriors = radon_2$posteriors) {
  agg <- dplyr::full_join(
    mcv_vals,
    axe_vals,
    by = c()
  ) %>%
    dplyr::full_join(posteriors) %>%
    dplyr::rename(elpd_cv_rst = rst_elpd) %>%
    dplyr::group_by(iter, perc, n_clusters, model, n_tot)


  agg
}



time_radon_simul <- function(info = radon_2) {
  dplyr::full_join(info$cv_yhats, info$posteriors) %>%
    dplyr::full_join(info$axe_yhats) %>%
    dplyr::group_by(iter, perc, n_clusters, model) %>%
    dplyr::select(dplyr::contains("time")) %>%
    unique() %>%
    dplyr::ungroup() %>%
    dplyr::summarise_at(dplyr::vars(dplyr::contains("time")), sum)
}

# figure_radon_simul <- function(results, info = radon_2) {
#
#     iter_plots <- purrr::map(info$n_clusters, function(nc) {
#         dt <-  results %>% dplyr::filter(n_clusters == nc)
#         lab <- unique(dt$nc_lab)
#         ggplot2::ggplot(
#             dt,
#             ggplot2::aes(x = cv_RMSE, y = axe_RMSE, color = factor(model))
#         ) +
#             ggplot2::geom_point(alpha = 0.4) +
#             ggplot2::facet_wrap(~ perc, scales = 'free', nrow = 1) +
#             ggplot2::geom_abline(slope = 1, intercept = 0, color = 'darkgray') +
#             ggplot2::geom_abline(slope = c(0.9, 1.1), intercept = 0,
#                                  color = 'gray', linetype = 'dashed') +
#             ggplot2::labs(x =  NULL, y = "", color = "Model") +
#             ggplot2::theme_bw() +
#             ggplot2::theme(
#                 legend.position = "bottom"
#             ) +
#             ggplot2::scale_color_brewer(type = "qual", palette = "Dark2") +
#             ggplot2::scale_x_continuous(breaks = scales::breaks_extended(n=4)) +
#             ggplot2::scale_y_continuous(breaks = scales::breaks_extended(n=4))
#
#
#
#     })
#
#     legend <- cowplot::get_legend(iter_plots[[1]])
#     iter_plots <- purrr::map(
#         iter_plots,
#         ~. + ggplot2::theme(legend.position = 'none')
#     )
#
#     ylab <-  cowplot::ggdraw() +
#         cowplot::draw_label(
#             expression('RMSE'['AXE']),
#             fontface = 'bold',
#             angle = 90
#         ) +
#         ggplot2::theme_bw() +
#         ggplot2::theme(panel.border = ggplot2::element_blank())
#     xlab <-  cowplot::ggdraw() +
#         cowplot::draw_label(
#             expression('RMSE'['MCV']),
#             fontface = 'bold'
#         ) +
#         ggplot2::theme_bw() +
#         ggplot2::theme(panel.border = ggplot2::element_blank())
#
#
#
#     cowplot::plot_grid(
#         cowplot::plot_grid(
#             ylab,
#             cowplot::plot_grid(
#                 plotlist = c(
#                     iter_plots,
#                     list(cowplot::plot_grid(NULL, legend, NULL, ncol = 1))),
#                 ncol = 2, align = "hv", axis = "tblr",
#                 labels = c(sprintf("J = %s", info$n_clusters), ""),
#                 hjust = -0.3, label_size = 10),
#             nrow = 1, rel_widths = c(0.02, 1)
#         ),
#         xlab, ncol = 1, rel_heights = c(1, 0.05)
#     )
# }

# unfortunately this is hard-coded to assume number of clusters is 3, 4, 6, 9, 12
# to change, adjust the the tidyselect::ends_with("#") lines to have the clusters
# of choice
# table_radon_simul <- function(results, latex = T) {
#     results %>%
#         dplyr::mutate_at(dplyr::vars(dplyr::starts_with("RMSE")), ~./cv_RMSE) %>%
#         dplyr::group_by(model, n_clusters) %>%
#         dplyr::summarise_at(
#             dplyr::vars(tidyselect::starts_with("RMSE")),
#             list(mean, sd)
#         ) %>%
#         dplyr::mutate_at(
#             dplyr::vars(tidyselect::starts_with("RMSE")),
#             round, digits = 2
#         ) %>%
#         tidyr::pivot_longer(
#             cols = dplyr::starts_with("RMSE"),
#             names_to = "method", values_to = "val"
#             ) %>%
#         dplyr::mutate(
#             statistic = dplyr::case_when(
#                 stringr::str_detect(method, "_fn1") ~ "ratio_avg",
#                 T ~ "ratio_sd"
#             )
#         )
#         tidyr::pivot_wider(id_cols = model, names_from = n_clusters,
#                            values_from = c(ratio_ave, ratio_sd)) %>%
#         dplyr::select(model,
#                       tidyselect::ends_with("3"),
#                       tidyselect::ends_with("4"),
#                       tidyselect::ends_with("6"),
#                       tidyselect::ends_with("9"),
#                       tidyselect::ends_with("12")
#         )
#
#
#
#     if (latex) {
#         # print the table in LaTeX
#         knitr::kable(tbl_ratios, "latex")
#     } else {
#         tbl_ratios
#     }
# }



#' @describeIn results_eight
#' @export
results_lol <- function(info = lol, posteriors = lol$posteriors,
                        axe = lol$axe_yhats, mcv_vals = lol$cv_yhats) {
  mcv <- do.call("rbind", mcv_vals)
  make_results_df(posteriors, axe, mcv) %>%
    dplyr::full_join(info$data %>%
      dplyr::mutate(idx = 1:dplyr::n()) %>%
      dplyr::select(idx, kills) %>%
      dplyr::rename(y = kills))
}


#' @describeIn results_eight
#' @export
results_slc <- function(info = slc, posteriors = slc$posteriors,
                        axe = slc$axe_yhats, mcv_vals = slc$cv_yhats) {
  make_results_df(posteriors, axe, mcv_vals) %>%
    dplyr::full_join(
      data.frame(
        y = info$data$y,
        idx = 1:length(info$data$y)
      )
    )
}



time_slc <- function(info = slc) {
  dplyr::full_join(info$cv_yhats, info$posteriors$df) %>%
    dplyr::group_by(loop) %>%
    dplyr::select(dplyr::contains("time")) %>%
    unique() %>%
    dplyr::ungroup() %>%
    dplyr::select(-loop) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::mutate(time_axe = info$axe_yhats$time)
}

make_results_df <- function(posteriors, axe, mcv_vals) {
  res <- dplyr::full_join(
    posteriors$df,
    axe$df,
    by = c()
  ) %>%
    dplyr::full_join(mcv_vals)

  # add AXE time
  J <- length(unique(res$loop))
  res %>%
    dplyr::group_by(loop) %>%
    dplyr::mutate(time_axe = axe$time / J)
}

#' @describeIn results_eight
#' @export
results_air <- function(info = air, posteriors = air$posteriors,
                        axe = air$axe_yhats, mcv_vals = air$cv_yhats) {
  make_results_df(posteriors, axe, mcv_vals) %>%
    dplyr::full_join(
      data.frame(
        y = info$data$df$observed,
        idx = 1:nrow(info$data$df)
      ),
      by = "idx"
    )
}





compare_time <- function(resdf) {
  resdf %>%
    dplyr::group_by(loop) %>%
    dplyr::select(dplyr::contains("time")) %>%
    unique() %>%
    dplyr::summarise_at(-loop, sum)
}

library(tidyverse)

# Loads into the nevironment all models and output in 'data' folder
for (obj in list.files("data")) {
  load(file.path("data", obj))
}
# Loads into the environment all functions in 'R' folder
for (file in list.files("R")) {
  source(file.path("R", file))
}

# convenient vector of data set names
datasets <- c("Eight schools", "Radon", "Radon subsets", "ESP", "SLC", "SRD")
methods <-  c("AXE", "GHOST", "VEHTARI", "iIS", "IJ-C",
              "NS-C", "IJ-A",  "NS-A", "PSIS-LOO")

base_colors <- c("gray30", RColorBrewer::brewer.pal(6, "Dark2"))
names(base_colors) <- c("AXE", "iIS", "GHOST", "IJ", "NS", "VEHTARI", "PSIS")
match_methods <- match(stringr::str_split(methods, "-", simplify = T)[, 1], names(base_colors))
method_colors <- base_colors[match_methods]
names(method_colors) <- methods
order <- c("AXE", "GHOST", "VEHTARI", "iIS", "IJ-C", "IJ-A", "NS-C", "NS-A")


# compile data frame with point-by-point, method-specific, and model-specific
# results in one data frame
acv_results <- function(subfolder) {
  acv_folder <- sprintf("python/results/%s", subfolder)
  fnms <- fs::dir_ls(acv_folder, type = 'file', regexp = "_yhats.csv")
  methods <- fnms |>
    basename() |> str_split("_") |> map_chr(function(x) pluck(x, 1)) |>
    str_replace("-c", "-a") |>
    map_chr(function(nm) {
      if (nm %in% c("IJ", "NS")) {
        nm <- paste(nm, "-c", sep = "")
      }
      nm
    })
  names(fnms) <- methods

  dfs <- map(methods, function(method) {
    suffix <-str_to_lower(method) |>  str_replace_all("-", "_")
    key <- list("yhat", "loop_time")
    names(key) <- c(sprintf("yhat_%s", suffix), sprintf("time_%s", suffix))
    read_csv(fnms[[method]]) |>
      rename(idx = i) |>
      select(-loop) |>
      mutate(yhat = as.numeric(yhat)) |>
      rename(!!!key)
  })
  df <- dfs[[1]]
  for (i in 2:length(dfs))
    df <- full_join(df, dfs[[i]])
  if (subfolder == "r2") {
    df <- mutate(df, perc = round(perc, digits = 1), idx = idx- 1)
  }

  df
}


eight_results <- results_eight() |>
  mutate(data_scale = round(data_scale, digits =1)) |>
  full_join(acv_results("eight_schools") |>
              mutate(data_scale = round(data_scale, digits = 1))) |>
  select(-contains("_iis_"), contains("time")) |>
  rename(yhat_psis_loo = yhat_psis_c) |>
  mutate(n = 1) |>
  group_by(loop, data_scale, n)


r1_results <- results_radon_full() |>
  full_join(acv_results("r1"), by = c("idx", "model")) |>
  group_by(model, loop) |>
  mutate(n = n()) |>
  group_by(model, loop, n) |>
  select(-contains("_iis_"), -yhat_axe2, contains("time")) # !!! Key step
r2_results <- results_radon_simul() |>
  mutate(perc = round(perc, digits = 1)) |>
  full_join(acv_results("r2")) |>
  rename(yhat_cv = cv_yhat) |>
  mutate(loop = interaction(model, perc, n_clusters, iter), n= 23) |>
  group_by(model, perc, n_clusters, iter, loop, n) |>
  select(-yhat_post_axe, -yhat_post, -contains("_iis_"), contains("time"))
lol_results <- results_lol() |>
  select(-contains("_iis_"), contains("time")) |>
  full_join(acv_results("lol")) |>
  group_by(loop) |>
  mutate(n = n()) |>
  group_by(loop, n) |>
  select(-yhat_post, -contains("ghst")) # |>
slc_results <- results_slc() |>
  select(-contains("_iis_"), -yhat_post, contains("time"), -contains("ghst")) |>
  full_join(acv_results("slc")) |>
  group_by(loop) |>
  mutate(n = n()) |>
  group_by(loop, n)
srd_results <- results_air() |>
  select(-yhat_post) |>
  select(-contains("_iis_"), contains("time"), -contains("ghst")) |>
  full_join(acv_results("air")) |>
  group_by(loop) |>
  mutate(n = n()) |>
  group_by(loop, n)

# group results df in a list
results_list <- list(eight_results, r1_results, r2_results,
                     lol_results, slc_results, srd_results)

# Create table 4 comparing results -----
# Calculates RMSE for each loop and combines all results outputs into one dataframe
rmse_all <- dplyr::bind_rows(
  purrr::map2(
    as.list(datasets),
    results_list,
    function(dname, dset) {
      if (!is.numeric(dset$loop)){
        dset$loop <- as.numeric(factor(dset$loop))
      }
      dset |>
        summarise_at(vars(starts_with("yhat_")), ~sqrt(mean((. - y)^2))) |>
        mutate(data = dname)
    })
) |>
  mutate(
    grouping = case_when(
      !is.na(data_scale) & data_scale < 2 ~ "alpha < 2",
      !is.na(data_scale) & data_scale >= 2 ~ "2 <= alpha",
      !is.na(model) ~ as.character(model),
      T ~ as.character(NA)
    )
  ) |>
  dplyr::select(data, n, loop, grouping,
                yhat_ns_c, yhat_ij_c, yhat_vt,
                yhat_ns_a, yhat_ij_a, yhat_psis_loo,
                yhat_psiis_fm, yhat_ghst_c,
                yhat_axe, yhat_cv, yhat_cv_sij)

vs_mean <- rmse_all |>
  # rename methods be consistent with other tables/figures
  mutate(ref = yhat_cv) |>
  mutate_at(vars(starts_with("yhat_")), ~abs(round(log(. / ref), digits = 3))) |>
  dplyr::select(-yhat_cv, -yhat_cv_sij, -starts_with("yhat_ij"), -starts_with("yhat_ns"))
vs_mode <- rmse_all |>
  # rename methods be consistent with other tables/figures
  select(data, loop, grouping, yhat_cv_sij, starts_with("yhat_ij"), starts_with("yhat_ns")) |>
  mutate(ref = yhat_cv_sij) |>
  mutate_at(vars(starts_with("yhat_")), ~round(log(./ref), digits = 3)) |>  dplyr::select(-yhat_cv_sij)


lrr_all <-
  full_join(vs_mean, vs_mode, by = c("data_scale", "data", "loop", "grouping")) |>
  pivot_longer(starts_with("yhat"), names_to = "statistic", values_to = "lrr") |>
  mutate(statistic = rename_methods(statistic) |>
           str_replace("YHAT-", "") |> factor(levels = order, ordered = TRUE),
         data = factor(data, levels = datasets, ordered = TRUE)) |>
  # calculate lrr mean, median, and sd
arrange(data, statistic) |>
    group_by(data, grouping, statistic) |>
  summarise(mean_lrr = mean(abs(lrr)),
            sd_lrr = sd(abs(lrr)),
            perc_under_25 = sum(abs(lrr) <= 0.25)/n()
            # med_lrr = median(lrr)
            ) |>
    # re-shape df and re-order columns for paper
  mutate_if(is.numeric, ~round(., digits = 2)) |>
  pivot_wider(id_cols = c(data, grouping),
              names_from = statistic,
              names_glue = "{statistic}_{.value}",
              values_from = c(mean_lrr, sd_lrr, perc_under_25)) |>
  tidyr::unite("data", data:grouping) |>
dplyr::select(data, contains("AXE"), contains("GHOST"), contains("VEHTARI"),
         contains("iIS-A"), contains("iIS"), contains("IJ"), contains("NS"))


dplyr::select(lrr_all, data,
               contains("mean"), contains("sd")) |>
        dplyr::select(data, contains("AXE"), contains("GHOST"),
                      contains("VEHTARI"), contains("iIS"),
                      contains("IJ-A"), contains("NS-A"))
# LaTeX output for LRR mean and standard deviations table
knitr::kable(dplyr::select(lrr_all, data,
                           contains("mean"), contains("sd")) |>
               dplyr::select(data, contains("AXE"), contains("GHOST"),
                             contains("VEHTARI"), contains("iIS"),
                             contains("IJ-A"), contains("NS-A")), "latex")


# Creates Figure 2 -------
# Combine yhat dataframes into one dataframe
yhats_all <- dplyr::bind_rows(
  purrr::map2(
    results_list,
    as.list(datasets),
    function(df, dname) {
      if (!('yhat_ij' %in% unique(df$method))) {
        df$yhat_ij <- NA
      }
      df |>
        dplyr::mutate(data = dname) |>
        dplyr::ungroup() |>
        dplyr::select(data, yhat_cv, yhat_axe, contains("yhat_ghst"), yhat_ij, contains("yhat_vt"), yhat_psiis_fm)
    }
  )
) |>
  tidyr::pivot_longer(
    cols = c(yhat_axe, contains("yhat_ghst"), yhat_psiis_fm, yhat_ij, yhat_vt),
    names_to = "method", values_to = "yhat"
  ) |>
  dplyr::mutate(
    method = stringr::str_replace(method, "yhat_", "") |>
      rename_methods(),
    data = factor(data, levels = datasets, ordered = T)
  )


# df with min/max values, to keep axes on the same scale for Figure 2
method_subset <-  c("AXE", "GHOST", "VEHTARI", "iIS")
yhats_minmax <- yhats_all |>
    dplyr::filter(method %in% method_subset[-3]) |>
    dplyr::group_by(data) |>
  dplyr::summarise(min = min(yhat), max = max(yhat), yhat_cv = yhat_cv[1]) |>
  tidyr::pivot_longer(cols = c(min, max), names_to = "range", values_to = "yhat") |>
  filter(!is.na(yhat))
# row A
ptbypt_axe <-
  yhats_all |>
  filter(method == "AXE") |>
  ggplot(aes(x = yhat_cv, y = yhat), color = method_colors["AXE"]) +
  theme_bw() +
  geom_point() +
  lemon::facet_rep_wrap(~data, scales = "free", nrow = 1) +
  labs(y = "AXE approximation", x = "MCV estimate") +
  geom_abline(slope = 1, intercept = 0) +
  theme(
    # text = element_text(family = "CMU Serif"),
    strip.background = element_rect(size = 0, fill = "white"),
    strip.text = element_text(hjust = 0, size = 12)
  ) +
  # invisible points at each data set's min/max to keep y-axis consistent
  geom_point(data = yhats_minmax, alpha = 0.01)



# row B
shapes <-  c(1, 2, 5, 0)
names(shapes) <- method_subset
p <- yhats_all |>
    dplyr::filter(method %in% method_subset[-3]) |>
    ggplot(aes(x = yhat_cv, y = yhat)) +
    geom_point(aes(color = method, shape = method)) +
    scale_color_manual(values = method_colors[method_subset[-3]]) +
    scale_shape_manual(values = shapes[-3]) +
    labs(
        y = "LCO Approximation", x = "MCV estimate",
        color = NULL, shape = NULL
    ) +
    lemon::facet_rep_wrap(~data, scales = "free", ncol = 6) +
    geom_abline(slope = 1, intercept = 0) +
    theme_bw() +
    theme(
        # text = element_text(family = "CMU Serif"),
        legend.key.size = unit(0.3, "lines"),
        legend.background = element_rect(
            fill = alpha("white", 0)
        ),
        strip.background = element_rect(size = 0, fill = "white"),
        strip.text = element_text(hjust = 0, size = 12)
    ) +
    geom_point(data = yhats_minmax, alpha = 0)

# combine rows A and B
p_both <- cowplot::plot_grid(
  ptbypt_axe,
  lemon::reposition_legend(
    p,
    panel = "panel-1-1", position = "top left", offset = c(0, -0.1)
  ),
  ncol = 1, labels = sprintf("%s)", LETTERS[1:2])
)

# saves Figure 2
ggsave("output/p_both.png", width = 10, height = 4.5, plot = p_both)



# Re-create figure 2, but including Vehtari ------

ptbypt_axe2 <-
  yhats_all |>
  filter(method == "AXE") |>
  ggplot(aes(x = yhat_cv, y = yhat), color = method_colors["AXE"]) +
  theme_bw() +
  geom_point() +
  lemon::facet_rep_wrap(~data, scales = "free", nrow = 1) +
  labs(y = "AXE approximation", x = "MCV estimate") +
  geom_abline(slope = 1, intercept = 0) +
  theme(
    # text = element_text(family = "CMU Serif"),
    strip.background = element_rect(size = 0, fill = "white"),
    strip.text = element_text(hjust = 0, size = 12)
  ) +
  # invisible points at each data set's min/max to keep y-axis consistent
  geom_point(data = yhats_minmax, alpha = 0)

p2 <- yhats_all |>
  dplyr::filter(method %in% method_subset) |>
  ggplot(aes(x = yhat_cv, y = yhat)) +
  geom_point(aes(color = method, shape = method)) +
  scale_color_manual(values = method_colors[method_subset]) +
    scale_shape_manual(values = shapes) +
  labs(
    y = "LCO Approximation", x = "MCV estimate",
    color = NULL, shape = NULL
  ) +
  lemon::facet_rep_wrap(~data, scales = "free", nrow = 1) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  theme(
    # text = element_text(family = "CMU Serif"),
    legend.key.size = unit(0.3, "lines"),
    legend.background = element_rect(
      fill = alpha("white", 0)
    ),
    strip.background = element_rect(size = 0, fill = "white"),
    strip.text = element_text(hjust = 0, size = 12)
  ) +
  geom_point(data = yhats_minmax, alpha = 0)

p_both2 <- cowplot::plot_grid(
  ptbypt_axe2,
  lemon::reposition_legend(
    p2,
    panel = "panel-1-1", position = "top left", offset = c(0, -0.1)
  ),
  ncol = 1, labels = sprintf("%s)", LETTERS[1:2])#,
  # label_fontfamily = "CMU Serif"
)

# saves Figure 2, including vehtari
ggsave("output/p_both2.png", width = 10, height = 4.5, plot = p_both2)


# radon subsets, limited to scenarios with at least 9 clusters
yhats_radon2 |>
  filter(n_clusters > 9) |>
  ggplot(aes(x = method, y = dif, fill = sprintf("Model %s", model))) +
  geom_boxplot() +
  theme_bw() +
  theme(
    legend.key.size = unit(0.75, "lines"),
    legend.background = element_rect(
      fill = alpha("white", 0)
    ),
    text = element_text(family = "CMU Serif")
  ) +
  labs(x = NULL, fill = NULL, y = NULL, title = "C) Radon subsets") +
  scale_fill_brewer(palette = "Set2")
# Not too different



# Times for all ------


# Randomly select 6 MCV folds
set.seed(39872)
mcv6 <- dplyr::bind_rows(
  purrr::map2(
    as.list(datasets[-3]),
    list(eight_results, r1_results, lol_results, slc_results, srd_results),
    function(dname, dset) {
      n_loops <- ifelse(dname == "Eight schools", 4, 6)
      sample_loops <- sample(unique(dset$loop), n_loops)
      if (dname == "Radon") {
        dset <- group_by(dset, model, loop)
      } else if (dname == "Eight schools") {
        message("Eight schools")
        dset <- group_by(dset, data_scale, loop)
      } else {
        dset <- group_by(dset, loop)
      }
      subset <- filter(dset, loop %in% sample_loops) |>
        mutate(time_cv = as.numeric(time_cv))
      subset |>
        summarise_at(vars(yhat_cv, yhat_axe), ~sqrt(mean((. - y)^2))) |>
        full_join(dplyr::select(subset, time_cv, loop) |> unique()) |>
        mutate(data = dname, loop = as.numeric(factor(loop)))
    })
) |> ungroup() |>
  mutate_at(vars(yhat_axe), ~log(./yhat_cv)) |>
 dplyr::select(-yhat_cv)


mcv6 |>
  group_by(data) |>
  summarise(mean_lrr = mean(abs(yhat_axe)), sd_lrr = sd(yhat_axe), time_cv = sum(time_cv))


# separate comparison for eight schools
subset_eight <- mcv6 |> filter(data == "Eight schools") |>
  group_by(data_scale) |>
  summarise(mean_lrr_subset = mean(abs(yhat_axe)), time_cv = sum(time_cv),
            sd_lrr = sd((yhat_axe)))


lrr_eight <- eight_results |>
  summarise_at(vars(starts_with("yhat_")), ~sqrt(mean((. - y)^2))) |>
  mutate(data = "Eight Schools") |>
  mutate(lrr_axe = log(yhat_axe/yhat_cv)) |>
  select(loop, data_scale, lrr_axe) |>
  group_by(data_scale) |> summarise(mean_lrr_axe = mean(abs(lrr_axe)))

full_join(
  subset_eight,
  lrr_eight
) |>
  mutate(high_diagn = abs(mean_lrr_subset) >= 0.25 | sd_lrr >= 0.25,
         high_axe = abs(mean_lrr_axe) >= 0.15) |>
  ungroup() |> select(high_diagn, high_axe) |> table()




# LRR figure 1 -----
method_subset <-  c("AXE", "GHOST",  "iIS")
dlabels <- sprintf("%s) %s", LETTERS[1:6], datasets)
lrr_all_vs_perc <- bind_rows(
  mutate(vs_mean, reference = 'mean'),
  mutate(vs_mode, reference = 'mode')
) |>
  pivot_longer(starts_with("yhat"), names_to = "method", values_to = "lrr") |>
  mutate(
    method = rename_methods(method) |>
           str_replace("YHAT-", ""),
    round_lrr = round(abs(lrr), digits = 2)
  ) |>
  group_by(data, method, round_lrr, reference) |>
  summarise(n = n())  |>
  group_by(data, method, reference) |>
  arrange(round_lrr) |>
  mutate(n_under = cumsum(n),
         perc_under = n_under/sum(n)) |>
  arrange(data, method, round_lrr) |>
  mutate(
    method = rename_methods(method),
    data = factor(data, levels = datasets, ordered = T) |>
      fct_relabel(~dlabels[match(., datasets)])
  )


max_lrr_all <- log(2)
j <- length(method_subset) - 1
lwd <- c(1.5, rep(0.75, j))
names(lwd) <- method_subset


p_lrr <- lrr_all_vs_perc |>
    filter(method %in% method_subset) |>
  filter(round_lrr <= max_lrr_all, !is.na(method)) |>
    mutate(method = factor(method, levels = method_subset, ordered = TRUE)) |>
  ggplot(aes(x = round_lrr, y = perc_under, color = method,
             group = interaction(method, reference))) +
  geom_line(aes(size = method, linetype = method)) +
  scale_size_manual(values = lwd) +
  scale_color_manual(values = method_colors[method_subset]) +
  facet_wrap(~data, nrow = 2) +
  labs(
    x = "|LRR|, truncated at log(2)",
       y =  expression(paste("Proportion of CV folds j with ", ~group("|", "LRR", "|")[j] <= group("|", "LRR", "|"))),
    color = NULL, size = NULL, linetype = NULL) +
  theme_bw() +
  theme(
    legend.background = element_rect(
      fill = alpha("transparent", 0)
    ),
    text = element_text(family = "CMU Serif"),
    strip.background = element_rect(size = 0, fill = "white"),
    strip.text = element_text(hjust = 0, size = 12)
  ) +
  guides(col = guide_legend(nrow = 3, byrow = TRUE))

p_lrr <- lemon::reposition_legend(p_lrr, panel = "panel-2-2", position = "bottom right")
ggsave("output/lrr_percentage.png", plot = p_lrr, width = 9.5, height = 5)


# LRR figure appendix -----
method_subset <- order
j <- length(method_subset) - 1
lwd <- c(1.5, rep(0.75, j))
method_lty <- map_chr(method_subset, ~ifelse(str_detect(., "-C"), "dashed", "solid"))
names(method_lty) <- method_subset
p_lrr_all <-
  lrr_all_vs_perc |>
  filter(round_lrr <= max_lrr_all, !is.na(method),
         method %in% method_subset) |>
  mutate(method = factor(method, levels = method_subset, ordered = TRUE)) |>
  ggplot(aes(x = round_lrr, y = perc_under, color = method)) +
  geom_line(aes(size = method, linetype = method)) +
  scale_linetype_manual(values = method_lty) +
  scale_size_manual(values = lwd) +
  scale_color_manual(values = method_colors[method_subset]) +
  facet_wrap(~data) +
  labs(
    x = "|LRR|, truncated at log(2)",
    y =  expression(paste("Proportion of CV folds j with ", ~group("|", "LRR", "|")[j] <= group("|", "LRR", "|"))),
    color = NULL, size = NULL, linetype = NULL) +
  theme_bw() +
  theme(
    # legend.position = 'bottom',
    legend.background = element_rect(
      fill = alpha("transparent", 0)
    ),
    strip.background = element_rect(size = 0, fill = "white"),
    strip.text = element_text(hjust = 0, size = 12)
  ) +
  guides(col = guide_legend(nrow = 3, byrow = TRUE))

# p_lrr_all <- lemon::reposition_legend(p_lrr_all, panel = "panel-3-2", position = "bottom right")
ggsave("output/lrr_percentage_all.png", plot = p_lrr_all, width = 9.5, height = 5)



# LRR AUCs ------
auc <- lrr_all_vs_perc |> filter(round_lrr <= max_lrr_all) |>
  group_by(data) |> mutate(max_lrr_all = max(round_lrr, na.rm = T)) |>
  group_by(data, method, max_lrr_all) |>
  arrange(round_lrr) |>
  mutate(
    width = c(0, round_lrr[-1] - round_lrr[-n()]),
    lhs = c(0, perc_under[-n()]),
    area = (lhs + perc_under)/2*width
  ) |>
  summarise(tot_area = sum(area), max_lrr = max(round_lrr)) |>
  mutate(auc = round((tot_area + (max_lrr_all - max_lrr))/max_lrr_all, digits = 2)) |>
    filter(method != "iIS-A")

auc_wide <- pivot_wider(auc, data, names_from = method, values_from = auc)
auc_wide <- auc_wide[, c(1:3, 10, 4:9)]
knitr::kable(auc_wide, "latex")

##
# Calculate times
##

times_all <- dplyr::bind_rows(
  purrr::map2(
    as.list(datasets),
    results_list,
    function(dname, dset) {
      if (!is.numeric(dset$loop)){
        dset$loop <- as.numeric(factor(dset$loop))
      }
      if ("axe_time" %in% colnames(dset))
        dset <- rename(dset, time_axe = axe_time)
      dset |>
        select(contains("time")) |>
        unique() |>
        mutate(data = dname) |>
        mutate_at(vars(contains("time")), as.numeric)
    })
)

times_df <- times_all |>
  group_by(data) |>
  summarise_at(vars(contains("time")), sum, na.rm = TRUE) |>
    mutate(data = factor(data, ordered = TRUE, levels = datasets)) |>
    arrange(data) |>
    mutate(data = dlabels)
names(times_df)[-1] <-  rename_methods(gsub("time_", "", names(times_df)[-1]))
times_df[, c(1, 4, 3, 2, 10, 6:9, 5)] |>
    mutate_if(is.numeric, round, digits = 1)



# A closer look at VT -----

ggplot(eight_results, aes(x = yhat_cv, y = yhat_vt, color = data_scale)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    labs(x = "Yhat (Manual CV)", y = "Yhat (Vehtari)") +
    theme_minimal()

ggplot(eight_results, aes(x = y, y = yhat_vt)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    theme_minimal() +
    labs(y = "Vehtari's method", x = "Y")



library(extrafont)
library(tidyverse)


for (file in list.files("R")) {
  source(file.path("R", file))
}

for (obj in list.files("data")) {
  load(file.path("data", obj))
}

eight_results <- results_eight() %>%
  select(-contains("_iis_")) %>%
  group_by(loop, data_scale)
r1_results <- results_radon_full() %>%
  group_by(model, loop) %>%
  select(-contains("_iis_"), -yhat_axe2) # !!! Key step
r2_results <- results_radon_simul() %>%
  rename(yhat_cv = cv_yhat) %>%
  mutate(loop = interaction(model, perc, n_clusters, iter)) %>%
  group_by(model, perc, n_clusters, iter, loop) %>%
  # rename(yhat_psiis_im = yhat_psiis_im.value) %>%
  select(-yhat_post_axe, -yhat_post, -contains("_iis_"))
lol_results <- results_lol() %>%
  select(-contains("_iis_")) %>%
  group_by(loop) %>%
  mutate(n = n()) %>%
  group_by(loop, n) %>%
  select(-yhat_post) # %>%
slc_results <- results_slc() %>%
  select(-contains("_iis_"), -yhat_post) %>%
  group_by(loop)
srd_results <- results_air() %>%
  select(-yhat_post) %>%
  select(-contains("_iis_")) %>%
  group_by(loop)



yhats_eight <- results_8 %>%
  df_compare_methods("yhat")
yhats_radon1 <- results_r1 %>%
  dplyr::mutate(n = dplyr::n()) %>%
  dplyr::group_by(model, loop, n) %>%
  df_compare_methods("yhat")
yhats_radon2 <- df_compare_methods(results_r2, "yhat")
# dplyr::mutate_at(dplyr::vars(dplyr::starts_with("yhat")), exp)
yhats_lol <- df_compare_methods(lol_results, "yhat")
yhats_slc <- df_compare_methods(results_slc, "yhat")
yhats_air <- df_compare_methods(results_air, "yhat")


# Figure 1
datasets <- c("Eight schools", "Radon", "Radon subsets", "ESP", "SLC", "SRD")
yhats_all <- dplyr::bind_rows(
  purrr::map2(
    list(results_8, results_r1, results_r2, lol_results, results_slc, results_air),
    as.list(datasets),
    function(df, dname) {
      df %>%
        dplyr::mutate(data = dname) %>%
        dplyr::ungroup() %>%
        dplyr::select(data, yhat_cv, yhat_axe, yhat_ghst_c, yhat_psiis_fm)
    }
  )
) %>%
  tidyr::pivot_longer(
    cols = c(yhat_axe, yhat_ghst_c, yhat_psiis_fm),
    names_to = "method", values_to = "yhat"
  ) %>%
  dplyr::mutate(
    method = stringr::str_replace(method, "yhat_", "") %>%
      rename_methods(),
    data = factor(data, levels = datasets, ordered = T)
  )


# df with min/max values
yhats_minmax <- yhats_all %>%
  dplyr::group_by(data) %>%
  dplyr::summarise(min = min(yhat), max = max(yhat), yhat_cv = yhat_cv[1]) %>%
  tidyr::pivot_longer(cols = c(min, max), names_to = "range", values_to = "yhat")

ptbypt_axe <-
  yhats_all %>%
  filter(method == "AXE") %>%
  ggplot(aes(x = yhat_cv, y = yhat)) +
  theme_bw() +
  geom_point() +
  lemon::facet_rep_wrap(~data, scales = "free", nrow = 1) +
  labs(y = "AXE approximation", x = "MCV estimate") +
  geom_abline(slope = 1, intercept = 0) +
  theme(
    text = element_text(family = "CMU Serif"),
    strip.background = element_rect(size = 0, fill = "white"),
    strip.text = element_text(hjust = 0, size = 12)
  ) +
  # invisible points at each data set's min/max to keep y-axis consistent
  geom_point(data = yhats_minmax, alpha = 0)



colors <- RColorBrewer::brewer.pal(6, "Dark2")[4:6]



p <- yhats_all %>%
  dplyr::filter(method != "iIS-A") %>%
  ggplot(aes(x = yhat_cv, y = yhat)) +
  geom_point(aes(color = method, shape = method)) +
  scale_color_manual(values = c("black", colors)) +
  labs(
    y = "LCO Approximation", x = "MCV estimate",
    color = NULL, shape = NULL
  ) +
  theme_bw() +
  lemon::facet_rep_wrap(~data, scales = "free", nrow = 1) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  theme(
    text = element_text(family = "CMU Serif"),
    legend.key.size = unit(0.3, "lines"),
    legend.background = element_rect(
      fill = alpha("white", 0)
    ),
    strip.background = element_rect(size = 0, fill = "white"),
    strip.text = element_text(hjust = 0, size = 12)
  ) +
  geom_point(data = yhats_minmax, alpha = 0)


p_both <- cowplot::plot_grid(
  ptbypt_axe,
  lemon::reposition_legend(
    p,
    panel = "panel-1-1", position = "top left", offset = c(0, -0.1)
  ),
  ncol = 1, labels = sprintf("%s)", LETTERS[1:2]),
  label_fontfamily = "CMU Serif"
)

ggsave("data-raw/p_all.png", width = 10, height = 6, plot = p_all)

ggsave("data-raw/p_both.png", width = 10, height = 4.5, plot = p_both)


# Figure 2
p_8 <- yhats_eight %>%
  mutate(
    alpha = case_when(
      data_scale < 1 ~ "low",
      data_scale < 2 ~ "mid",
      data_scale >= 2 ~ "high"
    ) %>%
      factor(levels = c("low", "mid", "high"), ordered = T) # ,
    # alpha = recode_factor(alpha,
    #                       low = "paste(phantom(12), alpha < 1)", mid = "1 <= {alpha < 2}",
    #                       high = "~~~alpha >= 2")
  ) %>%
  ggplot(aes(x = method, y = dif, fill = alpha)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    legend.key.size = unit(0.75, "lines"),
    legend.background = element_rect(
      fill = alpha("white", 0)
    ),
    text = element_text(family = "CMU Serif")
  ) +
  labs(x = NULL, fill = NULL, y = NULL, title = "A) Eight schools") +
  scale_fill_brewer(
    palette = "BuGn", type = "seq",
    # labels = scales::parse_format())+
    labels = c(
      expression(alpha < 1),
      expression(1 <= {
        alpha < 2
      }),
      expression(alpha >= 2),
      expression(123456)
    )
  ) +
  geom_hline(yintercept = 0)


p_8 <- lemon::reposition_legend(p_8, position = "bottom left")

p_r1 <- ggplot(
  yhats_radon1,
  aes(x = method, y = dif, fill = sprintf("Model %s", model))
) +
  geom_boxplot() +
  theme_bw() +
  theme(
    legend.key.size = unit(0.75, "lines"),
    legend.background = element_rect(
      fill = alpha("white", 0)
    ),
    text = element_text(family = "CMU Serif")
  ) +
  labs(x = NULL, fill = NULL, y = NULL, title = "B) Radon") +
  scale_fill_brewer(palette = "Set2") +
  geom_hline(yintercept = 0)
p_r1 <- lemon::reposition_legend(p_r1, position = "bottom left")


r2_cutoff <- data.frame(
  dif = -1,
  model = 3,
  method = c("GHOST", "iIS-C")
)
p_r2 <-
  yhats_radon2 %>%
  ggplot(aes(x = method, y = dif)) +
  geom_boxplot(aes(fill = sprintf("Model %s", model))) +
  ylim(-1, 1) +
  geom_point(
    aes(group = model),
    position = position_nudge(x = .25),
    data = r2_cutoff, shape = 3
  ) +
  theme_bw() +
  theme(
    legend.key.size = unit(0.75, "lines"),
    legend.background = element_rect(
      fill = alpha("white", 0)
    ),
    text = element_text(family = "CMU Serif")
  ) +
  labs(
    x = NULL, fill = NULL, y = NULL,
    title = "C) Radon subsets"
  ) +
  scale_fill_brewer(palette = "Set2") +
  geom_hline(yintercept = 0)
p_r2 <- lemon::reposition_legend(p_r2, position = "bottom left")

p_lol <- plot_compare_methods(yhats_lol, "yhat") +
  labs(y = NULL, title = "D) ESP")
p_slc <- plot_compare_methods(yhats_slc, "yhat") +
  labs(y = NULL, title = "E) SLC")
p_air <- plot_compare_methods(yhats_air, "yhat") +
  labs(y = NULL, subtitle = "F) SRD")

p_boxpl <- cowplot::plot_grid(p_8, p_r1, p_r2, p_lol, p_slc, p_air, nrow = 2)

ggsave("data-raw/compare_all.png", plot = p_boxpl, width = 10, height = 5)




# radon subsets, with more clusters
yhats_radon2 %>%
  filter(n_clusters > 9) %>%
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

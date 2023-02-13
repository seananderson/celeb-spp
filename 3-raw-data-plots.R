library(MASS)
library(dplyr)
library(ggplot2)
library(brms)

source("1-prep-data.R")
dir.create("figs", showWarnings = FALSE)

group_by(d, taxonomic_group, serial_number) %>%
  mutate(celeb_views = mean(celeb_average_daily_view, na.rm = TRUE)) %>%
  filter(celeb_views >= 1) %>%
  filter(species_average_daily_view > 0) %>%
  reframe(.diff = species_average_daily_view[celebrity == 1] / species_average_daily_view[celebrity == 0], celeb_views = mean(celeb_views, na.rm = TRUE)) %>%
  ggplot(aes(celeb_views, .diff)) +
  geom_point(alpha = 0.3) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~taxonomic_group, scales = "fixed") +
  ggsidekick::theme_sleek() +
  geom_hline(yintercept = 1, lty = 2) +
  geom_smooth(se = FALSE, method = "rlm", colour = "red", formula = y ~ x) +
  ylab("Ratio of average views\n(celebrity to non-celebrity species)") +
  xlab("Celebrity average page views")

ggsave("figs/raw-dat-rlm.png", width = 9, height = 5)

xx <- group_by(dpos, taxonomic_group, serial_number) %>%
  summarise(
    .diff = species_average_daily_view[celebrity == 1] / species_average_daily_view[celebrity == 0],
    celeb_views = mean(celeb_average_daily_view)
  ) %>%
  ungroup()
zz <- purrr::map_dfr(c(1, 10, 100, 1000), function(.x) {
  xx %>%
    filter(celeb_views > .x) %>%
    mutate(cutoff = .x)
})

zz %>%
  ggplot(aes(as.factor(cutoff), .diff)) +
  geom_point(position = position_jitter(height = 0, width = 0.1), alpha = 0.1) +
  geom_violin(draw_quantiles = c(0.5), alpha = 0.5, colour = "red") +
  scale_y_log10() +
  coord_cartesian(ylim = c(0.05, 20)) +
  ggsidekick::theme_sleek() +
  geom_hline(yintercept = 1, lty = 2) +
  ylab("Ratio of average views\n(celebrity to non-celebrity species)") +
  xlab("Celebrity daily views threshold") +
  facet_wrap(~taxonomic_group, scales = "fixed")

ggsave("figs/raw-dat-violin.png", width = 9, height = 5)

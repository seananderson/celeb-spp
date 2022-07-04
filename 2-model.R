library(dplyr)
library(ggplot2)
library(brms)

source("1-prep-data.R")
dir.create("figs", showWarnings = FALSE)

# get_prior(bf(
#   log(species_average_daily_view) ~ 1 + celebrity +
#     (1 + celebrity | taxonomic_group) + (1 | serial_number),
#   nu = 5
# ),
# data = dpos, family = student()
# )
#
# get_prior(bf(
#   species_average_daily_view ~ 1 + celebrity +
#     (1 + celebrity | taxonomic_group) + (1 | serial_number)
# ),
# data = dpos, family = negbinomial()
# )

fit_brms_mod <- function(dat, model_disp = FALSE, iter = 500L, chains = 4L,
                         family = c("Gamma", "Student", "NB2"), .nu = 5,
                         include_celeb_views = FALSE, independent_taxa = FALSE) {
  family <- match.arg(family)

  if (family == "Gamma") {
    .family <- Gamma(link = "log")
  }
  if (family == "Student") {
    .family <- student
  }
  if (family == "NB2") {
    .family <- brms::negbinomial()
  }

  priors <-
    prior(normal(0, 1), class = "b") +
    prior(normal(0, 5), class = "Intercept") +
    prior(student_t(3, 0, 2.5), class = "sd") +
    prior(lkj_corr_cholesky(1), class = "L")

  if (!independent_taxa) {
  if (!include_celeb_views) {
    if (family == "Gamma") {
      if (model_disp) {
        priors <- priors + prior(student_t(3, 0, 2.5), class = "b", dpar = "shape")
        f <- bf(
          species_average_daily_view ~
            1 + celebrity +
            (1 + celebrity | taxonomic_group) + (1 | serial_number),
          shape ~ 0 + taxonomic_group
        )
      } else {
        priors <- priors + prior(student_t(3, 0, 2.5), class = "shape")
        f <- bf(species_average_daily_view ~
          1 + celebrity +
          (1 + celebrity | taxonomic_group) + (1 | serial_number))
      }
    }

    if (family == "Student") {
      if (model_disp) {
        priors <- priors + prior(student_t(3, 0, 2.5), class = "b", dpar = "sigma")
          # prior(gamma(2, 0.1), class = "nu")
          # prior(gamma(4, 1), class = "nu")
        f <- bf(
          log(species_average_daily_view) ~
            1 + celebrity +
            (1 + celebrity | taxonomic_group) + (1 | serial_number),
          sigma ~ 0 + taxonomic_group, nu = .nu
        )
      } else {
        priors <- priors + prior(student_t(3, 0, 2.5), class = "sigma")
          # prior(gamma(2, 0.1), class = "nu")
          # prior(gamma(4, 1), class = "nu")
        f <- bf(log(species_average_daily_view) ~
          1 + celebrity +
          (1 + celebrity | taxonomic_group) + (1 | serial_number), nu = .nu)
      }
    }

    if (family == "NB2") {
      # priors <- priors + prior(gamma(0.01, 0.01), class = "shape")
      priors <- priors + prior(student_t(3, 0, 5), class = "shape")
      f <- bf(species_total_views ~
        1 + celebrity +
        (1 + celebrity | taxonomic_group) + (1 | serial_number))
    }
  } else {
    if (family == "NB2") {
      # priors <- priors + prior(gamma(0.01, 0.01), class = "shape")
      priors <- priors + prior(student_t(3, 0, 5), class = "shape")
      f <- bf(species_total_views ~
        1 + celebrity * log(celeb_average_daily_view) +
        (1 + celebrity * log(celeb_average_daily_view) | taxonomic_group) + (1 | serial_number))
    }
  }
  }

  if (independent_taxa) {
    priors <- prior(normal(0, 1), class = "b") +
      prior(normal(0, 5), class = "Intercept") +
      prior(student_t(3, 0, 2.5), class = "sd")
      # prior(lkj_corr_cholesky(1), class = "L")
      # prior(student_t(3, 0, 5), class = "shape")
    f <- bf(species_total_views ~
        1 + celebrity * taxonomic_group + (1 | serial_number), nu = .nu)
  }

  brm(
    formula = f,
    data = dat,
    family = .family,
    backend = "cmdstanr",
    iter = iter,
    chains = chains,
    cores = future::availableCores(logical = FALSE),
    prior = priors,
    control = list(adapt_delta = 0.98, max_treedepth = 12)
  )
}

ITER <- 2000L
CHAINS <- 4L

fit1_nb <- fit_brms_mod(d,
  iter = ITER, model_disp = FALSE,
  chains = CHAINS, family = "NB2"
)
fit1_10_nb <- fit_brms_mod(d_10,
  iter = ITER, model_disp = FALSE,
  chains = CHAINS, family = "NB2"
)
fit1_100_nb <- fit_brms_mod(d_100,
  iter = ITER, model_disp = FALSE,
  chains = CHAINS, family = "NB2"
)
fit1_1000_nb <- fit_brms_mod(d_1000,
  iter = ITER, model_disp = FALSE,
  chains = CHAINS, family = "NB2"
)

NU <- 4
fit1_st_nu4 <- fit_brms_mod(dpos,
  iter = ITER, model_disp = FALSE,
  chains = CHAINS,
  family = "Student", .nu = NU
)
fit1_10_st_nu4 <- fit_brms_mod(dpos_10,
  iter = ITER, model_disp = FALSE,
  chains = CHAINS,
  family = "Student", .nu = NU
)
fit1_100_st_nu4 <- fit_brms_mod(dpos_100,
  iter = ITER, model_disp = FALSE,
  chains = CHAINS,
  family = "Student", .nu = NU
)
fit1_1000_st_nu4 <- fit_brms_mod(dpos_1000,
  iter = ITER, model_disp = FALSE,
  chains = CHAINS,
  family = "Student", .nu = NU
)

# plotting ----------------------

m_nb2 <- list(fit1_nb, fit1_10_nb, fit1_100_nb, fit1_1000_nb)
m_st4 <- list(fit1_st_nu4, fit1_10_st_nu4, fit1_100_st_nu4, fit1_1000_st_nu4)

names(m_nb2) <- paste0(c(1, 10, 100, 1000))
names(m_st4) <- paste0(c(1, 10, 100, 1000))

dir.create("data-generated", showWarnings = FALSE)
saveRDS(m_nb2, file = "data-generated/nb2-models.rds")
saveRDS(m_st4, file = "data-generated/st4-models.rds")

m_nb2 <- readRDS("data-generated/nb2-models.rds")
m_st4 <- readRDS("data-generated/st4-models.rds")

col_lab <- "Celebrity\nviews\nthreshold"

# global coef plots ------------------------

p <- lapply(m_st4, brms::as_draws_df)
p2 <- purrr::map_dfr(p, function(.x) {
  data.frame(b_celebrity = .x$b_celebrity)
}, .id = "threshold")
ggplot(p2, aes(exp(b_celebrity), colour = threshold, fill = threshold)) +
  geom_vline(xintercept = 1, lty = 2, colour = "grey30") +
  geom_density(alpha = 0.1, colour = NA) +
  geom_density(fill = NA) +
  coord_cartesian(xlim = c(0.8, 2.5)) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_light() +
  scale_y_continuous(expand = expansion(mult = c(0, .02))) +
  labs(colour = col_lab, fill = col_lab, y = "Density", x = "Multiplicative celebrity effect\nacross taxonomic groups")
ggsave("figs/global-effect-st4.png", width = 7, height = 4)

saveRDS(p2, "data-generated/global-posterior-st4.rds")

p <- lapply(m_nb2, brms::as_draws_df)
p2 <- purrr::map_dfr(p, function(.x) {
  data.frame(b_celebrity = .x$b_celebrity)
}, .id = "threshold")
ggplot(p2, aes(exp(b_celebrity), colour = threshold, fill = threshold)) +
  geom_vline(xintercept = 1, lty = 2, colour = "grey30") +
  geom_density(alpha = 0.1, colour = NA) +
  geom_density(fill = NA) +
  coord_cartesian(xlim = c(0.8, 2.5)) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_light() +
  scale_y_continuous(expand = expansion(mult = c(0, .02))) +
  labs(colour = col_lab, fill = col_lab, y = "Density", x = "Multiplicative celebrity effect\nacross taxonomic groups")
ggsave("figs/global-effect-nb2.png", width = 7, height = 4)

saveRDS(p2, "data-generated/global-posterior-nb2.rds")

# posterior predictive plots ------------------

bayesplot::bayesplot_theme_set(theme_light())
g <- list()
XLAB <- xlab("Species page views")
YLAB <- ylab("Density")
g[[1]] <- pp_check(fit1_nb, ndraws = 25) + scale_x_log10() +
  ggtitle("Celebrity threshold: 1") + XLAB + YLAB
g[[2]] <- pp_check(fit1_10_nb, ndraws = 25) + scale_x_log10() +
  ggtitle("Celebrity threshold: 10") + XLAB + YLAB
g[[3]] <- pp_check(fit1_100_nb, ndraws = 25) + scale_x_log10() +
  ggtitle("Celebrity threshold: 100") + XLAB + YLAB
g[[4]] <- pp_check(fit1_1000_nb, ndraws = 25) + scale_x_log10() +
  ggtitle("Celebrity threshold: 1000") + XLAB + YLAB
cowplot::plot_grid(plotlist = g)
ggsave("figs/nb2-ppcheck.png", width = 10, height = 6)

XLIM <- xlim(-10, 10)
XLAB <- xlab("ln average daily species views")
g <- list()
g[[1]] <- pp_check(fit1_st_nu4, ndraws = 25) + XLIM + XLAB + YLAB +
  ggtitle("Celebrity threshold: 1")
g[[2]] <- pp_check(fit1_10_st_nu4, ndraws = 25) + XLIM + XLAB + YLAB+
  ggtitle("Celebrity threshold: 10")
g[[3]] <- pp_check(fit1_100_st_nu4, ndraws = 25) + XLIM + XLAB + YLAB+
  ggtitle("Celebrity threshold: 100")
g[[4]] <- pp_check(fit1_1000_st_nu4, ndraws = 25) + XLIM + XLAB + YLAB+
  ggtitle("Celebrity threshold: 1000")
cowplot::plot_grid(plotlist = g)
ggsave("figs/st4-ppcheck.png", width = 10, height = 6)

get_draws <- function(fit) {
  p1 <- tidybayes::spread_draws(fit, r_taxonomic_group[taxa, param])
  p2 <- tidybayes::spread_draws(fit, b_celebrity)
  filter(p1, param == "celebrity") |>
    left_join(p2)
}

# tidybayes::get_variables(fit1_10_nb_ind)
# plot_violins <- function(draws_df) {
#   ggplot(draws_df, aes(taxa, exp(b_celebrity + r_taxonomic_group))) +
#     geom_hline(yintercept = 1, lty = 2) +
#     geom_violin(fill = NA) +
#     coord_flip() +
#     xlab("") +
#     ylab("Effect") +
#     theme_light()
# }

# post2 <- readRDS("data-generated/nb2-models.rds")
# draws <- get_draws(post2[[2]])
# plot_violins(draws)

# library(DHARMa)
# post2 <- readRDS("data-generated/st4-models.rds")
# post2 <- readRDS("data-generated/nb2-models.rds")
# model <- post2[[4]]
# dat <- dpos_1000
# dat <- d_1000

# model_check <- createDHARMa(
#   simulatedResponse = t(posterior_predict(model, ndraws = 500L)),
#   observedResponse = dat$species_total_views,
#   fittedPredictedResponse = apply(t(posterior_epred(model)), 1L, mean),
#   integerResponse = TRUE
# )
# plot(model_check)
#
# model_check_pos <- createDHARMa(
#   simulatedResponse = t(posterior_predict(model, ndraws = 300L)),
#   observedResponse = log(dat$species_average_daily_view),
#   fittedPredictedResponse = apply(t(posterior_epred(model)), 1L, mean),
#   integerResponse = TRUE
# )
# plot(model_check_pos)

make_prob_table <- function(draws_df) {
  group_by(draws_df, taxa) |>
    mutate(theta = exp(b_celebrity + r_taxonomic_group)) |>
    summarise(
      prob_gt_one = mean(theta > 1),
      CI95_lwr = quantile(theta, 0.025), median = quantile(theta, 0.5),
      CI95_upr = quantile(theta, 0.975)
    )
}

tb <- purrr::map_dfr(m_nb2, function(.x) {
  get_draws(.x) %>% make_prob_table()
}, .id = "celeb")
saveRDS(tb, "data-generated/nb2-probs.rds")
knitr::kable(tb, digits = 2L)

tb <- purrr::map_dfr(m_st4, function(.x) {
  get_draws(.x) %>% make_prob_table()
}, .id = "celeb")
saveRDS(tb, "data-generated/st4-probs.rds")
knitr::kable(tb, digits = 2L)

dot_line_plot <- function(post) {
  post %>%
    # filter(celeb != "Celeb1000") %>%
    mutate(effect = exp(b_celebrity + r_taxonomic_group)) %>%
    group_by(taxa, celeb) %>%
    summarise(
      lwr = quantile(effect, probs = 0.05),
      upr = quantile(effect, probs = 0.95),
      lwr2 = quantile(effect, probs = 0.25),
      upr2 = quantile(effect, probs = 0.75),
      med = quantile(effect, probs = 0.5)
    ) %>%
    ggplot(aes(taxa, med, ymin = lwr, ymax = upr, colour = celeb)) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    geom_linerange(aes(ymin = lwr2, ymax = upr2), lwd = 1,
      position = position_dodge(width = 0.5)) +
    coord_flip() +
    xlab("") +
    ylab("Multiplicative effect") +
    theme_light() +
    scale_y_log10() +
    scale_colour_viridis_d(end = 0.9) +
    scale_fill_viridis_d(end = 0.9) +
    labs(colour = "Celebrity\nviews\nthreshold")
}

pp <- purrr::map_dfr(m_nb2, function(.x) {
  get_draws(.x)
}, .id = "celeb")
dot_line_plot(pp)
ggsave("figs/dot-line-1-1000-nb2.png", width = 5, height = 5)

pp <- purrr::map_dfr(m_st4, function(.x) {
  get_draws(.x)
}, .id = "celeb")
dot_line_plot(pp)
ggsave("figs/dot-line-1-1000-st4.png", width = 5, height = 5)

# independent models ------------------------------------------------

tax <- c("Amphibian", "Bird", "Fish", "Invertebrate", "Mammal", "Reptile")
thresholds <- c(1, 10, 100, 1000)
all <- expand.grid(tax = tax, thresholds = thresholds)

fit_ml_models <- function(.x, .y, student = TRUE, df = 4) {
  print(paste(.x, .y))
  if (.y == 1) dat <- dpos
  if (.y == 10) dat <- dpos_10
  if (.y == 100) dat <- dpos_100
  if (.y == 1000) dat <- dpos_1000

  dd <- dplyr::filter(dat, taxonomic_group == .x)
  dd$serial_number <- as.factor(dd$serial_number)
  fam <- if (student) sdmTMB::student(df = df) else gaussian()
  m <- sdmTMB::sdmTMB(log(species_average_daily_view) ~
      celebrity + (1 | serial_number),
    data = dd, family = fam, spatial = "off", reml = TRUE,
    silent = FALSE
  )
  cis <- sdmTMB::tidy(m, conf.int = TRUE)[2L,,drop=FALSE]
  cis50 <- sdmTMB::tidy(m, conf.int = TRUE,
    conf.level = 0.5)[2L,c("conf.low", "conf.high"),drop=FALSE]
  names(cis50) <- c("conf.low50", "conf.high50")
  data.frame(cis, cis50, taxa = .x, threshold = .y)
}

fits <- purrr::map2_dfr(all$tax, all$thresholds, function(.x, .y)
  fit_ml_models(.x, .y, TRUE, df = 4))

fits

fits %>% filter(!is.na(conf.low)) %>%
  ggplot(aes(taxa, exp(estimate),
    ymin = exp(conf.low), ymax = exp(conf.high),
    colour = as.factor(threshold))) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = exp(conf.low50), ymax = exp(conf.high50)), lwd = 1,
    position = position_dodge(width = 0.5)) +
  coord_flip() +
  xlab("") +
  ylab("Multiplicative effect") +
  theme_light() +
  scale_y_log10() +
  scale_colour_viridis_d(end = 0.9) +
  scale_fill_viridis_d(end = 0.9) +
  labs(colour = "Celebrity\nviews\nthreshold")

ggsave("figs/independent-sdmTMB-models-student-4.png", width = 5, height = 5)

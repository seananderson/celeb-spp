library(tidyverse)
library(brms)

source("1-prep-data.R")

# get_prior(bf(
#   species_average_daily_view ~ 1 + celebrity +
#     (1 + celebrity | taxonomic_group) + (1 | serial_number),
#   shape ~ 0 + taxonomic_group),
#   data = dpos, family = Gamma(link = "log"))
#
# get_prior(bf(
#   species_average_daily_view ~ 1 + celebrity +
#     (1 + celebrity | taxonomic_group) + (1 | serial_number)),
#   data = dpos, family = Gamma(link = "log"))

get_prior(bf(
  log(species_average_daily_view) ~ 1 + celebrity +
    (1 + celebrity | taxonomic_group) + (1 | serial_number), nu = 5),
  data = dpos, family = student())

get_prior(bf(
  species_average_daily_view ~ 1 + celebrity +
    (1 + celebrity | taxonomic_group) + (1 | serial_number)),
  data = dpos, family = negbinomial())



fit_brms_mod <- function(dat, model_disp = FALSE, iter = 500L, chains = 4L,
  family = c("Gamma", "Student", "NB2"), .nu = 5) {

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

  if (family == "Gamma") {
    if (model_disp) {
      priors <- priors + prior(student_t(3, 0, 2.5), class = "b", dpar = "shape")
      f <-  bf(species_average_daily_view ~
          1 + celebrity +
          (1 + celebrity | taxonomic_group) + (1 | serial_number),
        shape ~ 0 + taxonomic_group)
    } else {
      priors <- priors + prior(student_t(3, 0, 2.5), class = "shape")
      f <-  bf(species_average_daily_view ~
          1 + celebrity +
          (1 + celebrity | taxonomic_group) + (1 | serial_number))
    }
  }

  if (family == "Student") {
    if (model_disp) {
      priors <- priors + prior(student_t(3, 0, 2.5), class = "b", dpar = "sigma") +
        prior(gamma(2, 0.1), class = "nu")
          # prior(gamma(4, 1),   class = nu),
      f <-  bf(log(species_average_daily_view) ~
          1 + celebrity +
          (1 + celebrity | taxonomic_group) + (1 | serial_number),
        shape ~ 0 + taxonomic_group)
    } else {
      priors <- priors + prior(student_t(3, 0, 2.5), class = "sigma") +
        prior(gamma(2, 0.1), class = "nu")
      f <-  bf(log(species_average_daily_view) ~
          1 + celebrity +
          (1 + celebrity | taxonomic_group) + (1 | serial_number))
    }
  }


  if (family == "NB2") {
    priors <- priors + prior(gamma(0.01, 0.01), class = "shape")
    f <-  bf(species_total_views ~
        1 + celebrity +
        (1 + celebrity | taxonomic_group) + (1 | serial_number))
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
    control = list(adapt_delta = 0.90, max_treedepth = 12)
  )
}

fit1 <- fit_brms_mod(dpos, iter = 200, model_disp = FALSE)
fit2 <- fit_brms_mod(dpos, iter = 200, model_disp = TRUE)
# loo::loo(fit1, fit2)

fit_1000 <- fit_brms_mod(dpos_1000, iter = 300, model_disp = FALSE, chains = 1)
fit_1000_t <- fit_brms_mod(dpos_1000, iter = 200, model_disp = FALSE, chains = 1, family = "Student", .nu = 9)

fit1_t <- fit_brms_mod(dpos, iter = 250, model_disp = FALSE, chains = 1, .nu = 3, family = "Student")

fit1_nb <- fit_brms_mod(d, iter = 250, model_disp = FALSE, chains = 1, family = "NB2")

# plotting ----------------------

# shinystan::launch_shinystan(fit1)

pp_check(fit1_nb, ndraws = 100) + scale_x_log10()

pp_check(fit1_t, ndraws = 100)
pp_check(fit1_t, ndraws = 100) + coord_fixed(xlim = c(-10, 10))
pp_check(fit1_t, ndraws = 100) + xlim(-10, 10)
pp_check(fit1, ndraws = 100) + scale_x_log10()

pp_check(fit_1000_t, ndraws = 100)
pp_check(fit_1000_t, ndraws = 100) + xlim(-10, 10)
pp_check(fit_1000, ndraws = 100) + scale_x_log10()

p1t <- brms::posterior_predict(fit1_t, ndraws = 100)
dim(p1t[,])
mean(apply(p1t, 1, min))
mean(apply(p1t, 1, max))
range(log(dpos$species_average_daily_view))

p1 <- brms::posterior_predict(fit1, ndraws = 100)
mean(apply(log(p1), 1, min))
mean(apply(log(p1), 1, max))
range(log(dpos$species_average_daily_view))

p1nb <- brms::posterior_predict(fit1_nb, ndraws = 100)
mean(apply(p1nb, 1, min))
mean(apply(p1nb, 1, max))
range(dpos$species_average_daily_view)
# hist(p1nb[,2])
# hist(p1nb[,2])

get_draws <- function(fit) {
  p1 <- tidybayes::spread_draws(fit, r_taxonomic_group[taxa, param])
  p2 <- tidybayes::spread_draws(fit, b_celebrity)
  filter(p1, param == "celebrity") |>
    left_join(p2)
}

plot_violins <- function(draws_df) {
  ggplot(draws_df, aes(taxa, exp(b_celebrity + r_taxonomic_group))) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_violin(fill = NA) +
    coord_flip() + xlab("") + ylab("Rating") +
    theme_light()
}

make_prob_table <- function(draws_df) {
  group_by(draws_df, taxa) |>
    mutate(theta = exp(b_celebrity + r_taxonomic_group)) |>
    summarise(prob_gt_one = mean(theta > 1),
      CI95_lwr = quantile(theta, 0.025), median = quantile(theta, 0.5),
      CI95_upr = quantile(theta, 0.975)) |>
    knitr::kable(digits = 2L)
}

get_draws(fit_1000) %>% plot_violins()
get_draws(fit_1000_t) %>% plot_violins()

get_draws(fit1_nb) %>% plot_violins()
get_draws(fit1) %>% plot_violins()
get_draws(fit1_t) %>% plot_violins()

get_draws(fit1) %>% make_prob_table()
get_draws(fit1_t) %>% make_prob_table()
get_draws(fit_1000_t) %>% make_prob_table()

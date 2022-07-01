library(dplyr)
library(ggplot2)
library(brms)

source("1-prep-data.R")
dir.create("figs", showWarnings = FALSE)

get_prior(bf(
  log(species_average_daily_view) ~ 1 + celebrity +
    (1 + celebrity | taxonomic_group) + (1 | serial_number),
  nu = 5
),
data = dpos, family = student()
)

get_prior(bf(
  species_average_daily_view ~ 1 + celebrity +
    (1 + celebrity | taxonomic_group) + (1 | serial_number)
),
data = dpos, family = negbinomial()
)

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
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
}

# tibble(x = seq(from = 0, to = 60, by = .1)) %>%
#   expand(x, nesting(alpha = c(2, 4),
#     beta  = c(0.1, 1))) %>%
#   mutate(density = dgamma(x, alpha, beta),
#     group   = rep(letters[1:2], times = n() / 2)) %>%
#
#   # plot
#   ggplot(aes(x = x, ymin = 0, ymax = density,
#     group = group, fill = group)) +
#   geom_ribbon(size = 0, alpha = 3/4) +
#   scale_fill_viridis_d(option = "B", direction = -1,
#     begin = 1/3, end = 2/3) +
#   scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
#   scale_y_continuous(NULL, breaks = NULL) +
#   coord_cartesian(xlim = c(0, 50)) +
#   theme(panel.grid = element_blank(),
#     legend.position = "none")

ITER <- 1000L
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

# fit1_nb_ind <- fit_brms_mod(d,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS, family = "NB2", independent_taxa = TRUE
# )
# fit1_10_nb_ind <- fit_brms_mod(d_10,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS, family = "NB2", independent_taxa = TRUE
# )
# fit1_100_nb_ind <- fit_brms_mod(d_100,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS, family = "NB2", independent_taxa = TRUE
# )
# fit1_1000_nb_ind <- fit_brms_mod(d_1000,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS, family = "NB2", independent_taxa = TRUE
# )

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



# fit1_st_nu4_ind <- fit_brms_mod(dpos,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS, independent_taxa = TRUE,
#   family = "Student", .nu = NU
# )
# fit1_10_st_nu4_ind <- fit_brms_mod(dpos_10,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS,independent_taxa = TRUE,
#   family = "Student", .nu = NU
# )
# fit1_100_st_nu4_ind <- fit_brms_mod(dpos_100,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS,independent_taxa = TRUE,
#   family = "Student", .nu = NU
# )
# fit1_1000_st_nu4_ind <- fit_brms_mod(dpos_1000,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS,independent_taxa = TRUE,
#   family = "Student", .nu = NU
# )

# glmmTMB basics??

# library(glmmTMB)

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
  # if (.y == 10 && .x == "Mammal") browser()
  m <- sdmTMB::sdmTMB(log(species_average_daily_view) ~
      celebrity + (1 | serial_number),
    data = dd, family = fam, spatial = "off", reml = TRUE,
    silent = FALSE
  )
  cis <- sdmTMB::tidy(m, conf.int = TRUE)[2L,,drop=FALSE]
  cis50 <- sdmTMB::tidy(m, conf.int = TRUE, conf.level = 0.5)[2L,c("conf.low", "conf.high"),drop=FALSE]
  names(cis50) <- c("conf.low50", "conf.high50")
  data.frame(cis, cis50, taxa = .x, threshold = .y)
}

# fits3 <- purrr::map2_dfr(all$tax, all$thresholds, function(.x, .y)
#   fit_ml_models(.x, .y, TRUE, df = 3))

fits <- purrr::map2_dfr(all$tax, all$thresholds, function(.x, .y)
  fit_ml_models(.x, .y, TRUE, df = 4))

# fits7 <- purrr::map2_dfr(all$tax, all$thresholds, function(.x, .y)
#   fit_ml_models(.x, .y, TRUE, df = 7))

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

# NU <- 5
# fit1_st_nu5 <- fit_brms_mod(dpos,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS,
#   family = "Student", .nu = NU
# )
# fit1_10_st_nu5 <- fit_brms_mod(dpos_10,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS,
#   family = "Student", .nu = NU
# )
# fit1_100_st_nu5 <- fit_brms_mod(dpos_100,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS,
#   family = "Student", .nu = NU
# )
# fit1_1000_st_nu5 <- fit_brms_mod(dpos_1000,
#   iter = ITER, model_disp = FALSE,
#   chains = CHAINS,
#   family = "Student", .nu = NU
# )

# fit1_nb_views <- fit_brms_mod(d, iter = 400,
# model_disp = FALSE, chains = 5, family = "NB2", include_celeb_views = TRUE)

# plotting ----------------------

# shinystan::launch_shinystan(fit1_nb)
# shinystan::launch_shinystan(fit1_st_nu4)

m_nb2 <- list(fit1_nb, fit1_10_nb, fit1_100_nb, fit1_1000_nb)
m_st4 <- list(fit1_st_nu4, fit1_10_st_nu4, fit1_100_st_nu4, fit1_1000_st_nu4)
# m_st4_ind <- list(fit1_st_nu4_ind, fit1_10_st_nu4_ind, fit1_100_st_nu4_ind, fit1_1000_st_nu4_ind)
# m_st5 <- list(fit1_st_nu5, fit1_10_st_nu5, fit1_100_st_nu5, fit1_1000_st_nu5)

names(m_nb2) <- paste0(c(1, 10, 100, 1000))
names(m_st4) <- paste0(c(1, 10, 100, 1000))
# names(m_st5) <- paste0(c(1, 10, 100, 1000))

dir.create("data-generated", showWarnings = FALSE)
saveRDS(m_nb2, file = "data-generated/nb2-models.rds")
saveRDS(m_st4, file = "data-generated/st4-models.rds")

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

# pp_check(fit1_st_nu5) + XLIM
# pp_check(fit1_10_st_nu5) + XLIM
# pp_check(fit1_100_st_nu5) + XLIM
# pp_check(fit1_1000_st_nu5) + XLIM

# pp_check(fit_1000_t, ndraws = 100)
# pp_check(fit_1000_t, ndraws = 100) + xlim(-10, 10)
# pp_check(fit_100_t, ndraws = 100) + xlim(-10, 10)
# pp_check(fit_1000, ndraws = 100) + scale_x_log10()

# p1t <- brms::posterior_predict(fit1_100_st_nu5, ndraws = 100)
# # dim(p1t[, ])
# mean(apply(p1t, 1, min))
# mean(apply(p1t, 1, max))
# range(log(dpos$species_average_daily_view))
#
# p1 <- brms::posterior_predict(fit2, ndraws = 100)
# mean(apply(log(p1), 1, min))
# mean(apply(log(p1), 1, max))
# range(log(dpos$species_average_daily_view))

# p1nb <- brms::posterior_predict(fit1_100_nb, ndraws = 100)
# mean(apply(p1nb, 1, min))
# mean(apply(p1nb, 1, max))
# range(dpos$species_average_daily_view)
# hist(log(p1nb[2,] + 1))
# hist(log(d$species_total_views + 1))
# # hist(p1nb[,2])

# z <- tidybayes::get_variables(fit1_nb)
# z[grep("logceleb", z)]

get_draws <- function(fit) {
  p1 <- tidybayes::spread_draws(fit, r_taxonomic_group[taxa, param])
  p2 <- tidybayes::spread_draws(fit, b_celebrity)
  # p3 <- tidybayes::spread_draws(fit, b_logceleb_average_daily_view)
  # p3 <- tidybayes::spread_draws(fit, b_logceleb_average_daily_view)
  filter(p1, param == "celebrity") |>
    left_join(p2)
}

# get_draws_ind <- function(fit) {
#   p <- tidybayes::spread_draws(fit, `b_.*`, regex = TRUE)
#   bird <-    p$`b_celebrity:taxonomic_groupBird`
#   invert <-  p$`b_celebrity:taxonomic_groupInvertebrate`
#   mammal <-  p$`b_celebrity:taxonomic_groupMammal`
#   reptile <- p$`b_celebrity:taxonomic_groupReptile`
#   fish <-    p$`b_celebrity:taxonomic_groupFish`
#   amphibian <- rep(0, nrow(p))
#   post <- bind_rows(
#     data.frame(taxa = "Bird", r_taxonomic_group = bird, stringsAsFactors = FALSE),
#     data.frame(taxa = "Invertebrate", r_taxonomic_group = invert, stringsAsFactors = FALSE),
#     data.frame(taxa = "Mammal", r_taxonomic_group = mammal, stringsAsFactors = FALSE),
#     data.frame(taxa = "Reptile", r_taxonomic_group = reptile, stringsAsFactors = FALSE),
#     data.frame(taxa = "Fish", r_taxonomic_group = fish, stringsAsFactors = FALSE),
#     data.frame(taxa = "Amphibian", r_taxonomic_group = amphibian, stringsAsFactors = FALSE),
#   )
#   post$b_celebrity <- rep(p$b_celebrity, 6L)
#   post
# }

# dd <- get_draws_ind(fit1_100_st_nu4_ind)


# tidybayes::get_variables(fit1_10_nb_ind)

plot_violins <- function(draws_df) {
  ggplot(draws_df, aes(taxa, exp(b_celebrity + r_taxonomic_group))) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_violin(fill = NA) +
    coord_flip() +
    xlab("") +
    ylab("Effect") +
    theme_light()
}

post2 <- readRDS("data-generated/nb2-models.rds")
draws <- get_draws(post2[[2]])
plot_violins(draws)
# get_draws_ind(fit1_10_nb_ind) %>% plot_violins()

# pred <- brms::posterior_predict(fit1_10_nb_ind, ndraws = 1)

library(DHARMa)
post2 <- readRDS("data-generated/st4-models.rds")
post2 <- readRDS("data-generated/nb2-models.rds")
model <- post2[[4]]
dat <- dpos_1000
dat <- d_1000

model_check <- createDHARMa(
  simulatedResponse = t(posterior_predict(model, ndraws = 500L)),
  observedResponse = dat$species_total_views,
  fittedPredictedResponse = apply(t(posterior_epred(model)), 1L, mean),
  integerResponse = TRUE
)
plot(model_check)

model_check_pos <- createDHARMa(
  simulatedResponse = t(posterior_predict(model, ndraws = 300L)),
  observedResponse = log(dat$species_average_daily_view),
  fittedPredictedResponse = apply(t(posterior_epred(model)), 1L, mean),
  integerResponse = TRUE
)
plot(model_check_pos)

make_prob_table <- function(draws_df) {
  group_by(draws_df, taxa) |>
    mutate(theta = exp(b_celebrity + r_taxonomic_group)) |>
    summarise(
      prob_gt_one = mean(theta > 1),
      CI95_lwr = quantile(theta, 0.025), median = quantile(theta, 0.5),
      CI95_upr = quantile(theta, 0.975)
    )
}

# get_draws(fit_1000) %>% plot_violins()
# get_draws(fit_1000_t) %>% plot_violins()
# get_draws(fit_100_t) %>% plot_violins()

# pp <- purrr::map_dfr(m_nb2, function(.x) {
#   get_draws(.x)
# }, .id = "celeb")

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

# pp %>%
#   # filter(celeb != "Celeb1000") %>%
#   ggplot(aes(taxa, exp(b_celebrity + r_taxonomic_group), colour = celeb, fill = celeb)) +
#   geom_hline(yintercept = 1, lty = 2) +
#   geom_violin() +
#   coord_flip() +
#   xlab("") +
#   ylab("Effect") +
#   theme_light() +
#   scale_y_log10() +
#   scale_colour_viridis_d(end = 0.9) +
#   scale_fill_viridis_d(end = 0.9)

# viol_dat <- pp %>% filter(celeb != "Celeb1000") %>%
# mutate(effect = exp(b_celebrity + r_taxonomic_group))

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
    # geom_violin(data = viol_dat, mapping = aes(x = taxa, y = effect, colour = celeb, fill = celeb ), inherit.aes = FALSE) +
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

# pp <- purrr::map_dfr(m_st4_ind, function(.x) {
#   get_draws_ind(.x)
# }, .id = "celeb")
# dot_line_plot(pp)

# get_draws(fit1_nb) %>% plot_violins()
# get_draws(fit1_10_nb) %>% plot_violins()
# get_draws(fit1_100_nb) %>% plot_violins()
# get_draws(fit1) %>% plot_violins()
# get_draws(fit2) %>% plot_violins()
# get_draws(fit1_t) %>% plot_violins()
# get_draws(fit2_t) %>% plot_violins()

# get_draws(fit1) %>% make_prob_table()
# get_draws(fit1_t) %>% make_prob_table()
# get_draws(fit_1000_t) %>% make_prob_table()

# tg <- select(d, taxonomic_group) %>% distinct()
#
# cv <- exp(seq(log(min(d$celeb_average_daily_view)), log(max(d$celeb_average_daily_view)), length.out = 100))
#
# nd <- expand.grid(
#   taxonomic_group = tg$taxonomic_group,
#   celebrity = c(0, 1),
#   celeb_average_daily_view = cv
# )
# nd <- arrange(nd, celebrity, taxonomic_group, celeb_average_daily_view)
# nrow(nd)
#
# p <- brms::posterior_linpred(fit1_nb_views, newdata = nd, re_formula = ~ (1 + celebrity * log(celeb_average_daily_view) | taxonomic_group), transform = TRUE)
# p_post <- brms::posterior_predict(fit1_nb_views, newdata = nd, re_formula = ~ (1 + celebrity * log(celeb_average_daily_view) | taxonomic_group))
# dim(p)
#
# nd$est <- apply(p, 2, median)
# nd$lwr <- apply(p, 2, quantile, probs = 0.1)
# nd$upr <- apply(p, 2, quantile, probs = 0.9)
#
# ggplot(nd, aes(celeb_average_daily_view, est, ymin = lwr, ymax = upr, colour = taxonomic_group, fill = taxonomic_group)) +
#   geom_ribbon(alpha = 0.2) +
#   geom_line() +
#   facet_grid(taxonomic_group ~ celebrity, scales = "free") +
#   scale_x_log10() +
#   # scale_y_log10() +
#   ggsidekick::theme_sleek()

# .n <- filter(nd, celebrity == 0) %>% nrow()
# ratio <- p[, (.n + 1):nrow(nd)] / p[, 1:.n]
#
# nd2 <- filter(nd, celebrity == 1) # pick one
# nd2$est <- apply(ratio, 2, median)
# nd2$lwr <- apply(ratio, 2, quantile, probs = 0.05)
# nd2$upr <- apply(ratio, 2, quantile, probs = 0.95)
#
# nd2$lwr2 <- apply(ratio, 2, quantile, probs = 0.25)
# nd2$upr2 <- apply(ratio, 2, quantile, probs = 0.75)

# nd2 %>%
#   # filter(celeb_average_daily_view >= 1) %>%
#   ggplot(aes(celeb_average_daily_view / 1000, est, ymin = lwr, ymax = upr)) +
#   geom_ribbon(alpha = 0.2) +
#   geom_ribbon(mapping = aes(ymin = lwr2, ymax = upr2), alpha = 0.2) +
#   geom_line() +
#   facet_wrap(~taxonomic_group, scales = "free_y") +
#   # coord_cartesian(ylim = c(1, NA), expand = FALSE)
#   # scale_y_continuous(lim = c(0.8, NA)) +
#   geom_hline(yintercept = 1, lty = 2) +
#   # scale_x_log10() +
#   # scale_y_log10() +
#   scale_x_continuous(trans = "sqrt", breaks = c(1, 10, 25, seq(50, 200, 50))) +
#   ggsidekick::theme_sleek() +
#   ylab("Multiplicative effect on species page views") +
#   xlab("Celebrity page views (1000s)")
#
# fit1_nb


# discrete celeb views? --------------------

# d$celeb_over10 <- as.integer(d$celeb_average_daily_view > 10)
# d$celeb_over100 <- as.integer(d$celeb_average_daily_view > 100)
# d$celeb_over1000 <- as.integer(d$celeb_average_daily_view > 1000)
#
# priors <-
#   prior(normal(0, 2), class = "b") +
#   prior(normal(0, 5), class = "Intercept") +
#   prior(student_t(3, 0, 2.5), class = "sd") +
#   prior(lkj_corr_cholesky(1), class = "L")
# priors <- priors + prior(student_t(3, 0, 2.5), class = "shape")
#
# f <- bf(species_total_views ~
#   1 + celebrity + celeb_over10 * celebrity + celeb_over100 * celebrity + celeb_over1000 * celebrity +
#   (1 + celebrity + celeb_over10 * celebrity + celeb_over100 * celebrity + celeb_over1000 * celebrity | taxonomic_group) + (1 | serial_number))
#
# fit <- brm(
#   formula = f,
#   data = d,
#   family = negbinomial(),
#   backend = "cmdstanr",
#   iter = 300,
#   chains = 1,
#   cores = future::availableCores(logical = FALSE),
#   prior = priors,
#   control = list(adapt_delta = 0.90, max_treedepth = 12)
# )
#
# nd <- expand.grid(
#   taxonomic_group = tg$taxonomic_group,
#   celebrity = c(0, 1)
# )
# nd <- bind_rows(
#   mutate(nd, celeb_over10 = 0, celeb_over100 = 0, celeb_over1000 = 0),
#   mutate(nd, celeb_over10 = 1, celeb_over100 = 0, celeb_over1000 = 0),
#   mutate(nd, celeb_over10 = 1, celeb_over100 = 1, celeb_over1000 = 0),
#   mutate(nd, celeb_over10 = 1, celeb_over100 = 1, celeb_over1000 = 1)
# )
#
# nd <- arrange(
#   nd, celebrity, taxonomic_group,
#   celeb_over10,
#   celeb_over100,
#   celeb_over1000,
# )
# nrow(nd)
#
# p <- brms::posterior_linpred(fit,
#   newdata = nd, re_formula = ~ (1 + celebrity + celeb_over10 * celebrity + celeb_over100 * celebrity + celeb_over1000 * celebrity | taxonomic_group),
#   transform = TRUE
# )
# dim(p)
#
# nd$est <- apply(p, 2, median)
# nd$lwr <- apply(p, 2, quantile, probs = 0.1)
# nd$upr <- apply(p, 2, quantile, probs = 0.9)

# ggplot(nd, aes(celeb_average_daily_view, est, ymin = lwr, ymax = upr, colour = taxonomic_group, fill = taxonomic_group)) +
#   geom_ribbon(alpha = 0.2) +
#   geom_line() +
#   facet_grid(taxonomic_group~celebrity, scales = "free") +
#   scale_x_log10()
#
# .n <- filter(nd, celebrity == 0) %>% nrow()
# ratio <- p[, (.n + 1):nrow(nd)] / p[, 1:.n]
#
# nd2 <- filter(nd, celebrity == 1) # pick one
# nd2$est <- apply(ratio, 2, median)
# nd2$lwr <- apply(ratio, 2, quantile, probs = 0.05)
# nd2$upr <- apply(ratio, 2, quantile, probs = 0.95)
# nd2$lwr2 <- apply(ratio, 2, quantile, probs = 0.25)
# nd2$upr2 <- apply(ratio, 2, quantile, probs = 0.75)
#
# nd2 <- mutate(nd2, celeb_cat = case_when(
#   celeb_over10 == 0 & celeb_over100 == 0 & celeb_over1000 == 0 ~ "Celeb_0001",
#   celeb_over10 == 1 & celeb_over100 == 0 & celeb_over1000 == 0 ~ "Celeb_0010",
#   celeb_over10 == 1 & celeb_over100 == 1 & celeb_over1000 == 0 ~ "Celeb_0100",
#   celeb_over10 == 1 & celeb_over100 == 1 & celeb_over1000 == 1 ~ "Celeb_1000"
# ))
# ggplot(nd2, aes(taxonomic_group, est, ymin = lwr, ymax = upr, colour = celeb_cat)) +
#   geom_pointrange(position = position_dodge(width = 0.5)) +
#   coord_flip() +
#   scale_y_log10() +
#   scale_color_viridis_d(end = 0.9) +
#   theme_light()



# raw data plots? -------------


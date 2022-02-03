library(tidyverse)
library(glmmTMB)
library(DHARMa)

# 8000 strong dataset, half of which consists of species named after celebrity
# (defined as a person who has a Wikipedia page with an average of 1 visit a day
# for the last 6 years) and the other half of their closest phylogenetic
# relative not named after a celebrity. Then I looked at average Wikipedia views
# for all of those species.
# d <- readr::read_csv(file = "data/raw/Celeb_Species_Reorg.csv")
# d$Amphibian <- NULL
# d$Bird <- NULL
# d$Fish <- NULL
# d$Mammal <- NULL
# d$Invertebrate <- NULL
# d$Reptile <- NULL
# saveRDS(d, file = "data/raw/Celeb_Species_Reorg.rds")

d <- readRDS("data/raw/Celeb_Species_Reorg.rds")
glimpse(d)
# Celebrity
# Celeb_Total_Views
# Celeb_Average_Daily_View
# Species_Total_Views
# Species_Average_Daily_View

names(d) <- tolower(names(d))
d$taxonomic_group <- as.factor(d$taxonomic_group)
# d <- select(d, -amphibian, -bird, -fish, -invertebrate, -mammal, -reptile)

ggplot(d, aes(as.factor(celebrity), log10(species_average_daily_view + 1))) +
  geom_boxplot()

ggplot(d, aes(as.factor(celebrity), log10(species_average_daily_view + 1))) +
  geom_violin()

ggplot(d, aes(as.factor(celebrity), log10(species_average_daily_view + 1))) +
  geom_violin() +
  facet_wrap(vars(taxonomic_group))

table(d$taxonomic_group)

ggplot(d, aes(as.factor(celebrity), log10(species_average_daily_view + 1))) +
  geom_violin() +
  facet_wrap(vars(taxonomic_group), scales = "free_y")

ggplot(d, aes(as.factor(celebrity), log10(species_average_daily_view + 1))) +
  geom_boxplot() +
  facet_wrap(vars(taxonomic_group), scales = "free_y")

# # naive, not paired:
# m <- glmmTMB(
#   species_total_views ~ 1 + celebrity +
#     (1 + celebrity | taxonomic_group),
#   data = d, family = glmmTMB::nbinom2(), verbose = TRUE
# )
# summary(m)
# ranef(m)
# coef(m)

# all same time for views?
range(d$celeb_total_views / d$celeb_average_daily_view, na.rm = TRUE)
range(d$species_total_views / d$species_average_daily_view, na.rm = TRUE)
x <- d$species_total_views / d$species_average_daily_view
x <- x / 365.25
x <- x[is.finite(x)]
x <- x[!is.na(x)]
hist(x, breaks = 100, main = "(celeb_total_views/celeb_average_daily_view)/365")

table(d$serial_number) %>% unique()
table(d$taxonomic_group)

# paired:
m1 <- glmmTMB(
  species_total_views ~ 1 + celebrity +
    (1 + celebrity | taxonomic_group) + (1 | serial_number),
  data = d, family = glmmTMB::nbinom2(), verbose = TRUE
)
m1$sdr
summary(m1)
b <- coef(m1)
b$cond$taxonomic_group

# b <- broom.mixed::tidy(m1, effects = "ran_vals")
# filter(b, term == "celebrity") %>%
#   ggplot(aes(estimate,
#     y = level,
#     xmin = estimate - 2 * std.error, xmax = estimate + 2 * std.error
#   )) +
#   geom_pointrange() +
#   geom_vline(xintercept = 0, lty = 2)
#
# res1 <- DHARMa::simulateResiduals(m1, plot = TRUE)

# maybe with dispersion varying by taxa?
m1.1 <- glmmTMB(
  species_total_views ~ 1 + celebrity +
    (1 + celebrity | taxonomic_group) + (1 | serial_number),
  dispformula = ~ 0 + taxonomic_group,
  data = d, family = glmmTMB::nbinom2(), verbose = TRUE
)
m1.1$sdr
summary(m1.1)
res1.1 <- DHARMa::simulateResiduals(m1.1, plot = TRUE)
AIC(m1, m1.1)
b <- broom.mixed::tidy(m1.1, effects = "ran_vals")
filter(b, term == "celebrity") %>%
  mutate(estimate = estimate + fixef(m1.1)$cond[[2]]) %>%
  ggplot(aes(exp(estimate),
    y = level,
    xmin = exp(estimate - 2 * std.error), xmax = exp(estimate + 2 * std.error)
  )) +
  geom_pointrange() +
  geom_vline(xintercept = 1, lty = 2)

# s <- simulate(m1.1, nsim = 500)

m1.2 <- glmmTMB(
  species_average_daily_view ~ 1 + celebrity +
    (1 + celebrity | taxonomic_group) + (1 | serial_number),
  dispformula = ~ 0 + taxonomic_group,
  data = d, family = glmmTMB::tweedie(), verbose = TRUE
)
summary(m1.2)
res1.2 <- DHARMa::simulateResiduals(m1.2, plot = TRUE)

# m2 <- glmmTMB(
#   species_total_views ~ 1 + celebrity +
#     (1 + celebrity | taxonomic_group) + (1 | serial_number),
#   data = d, family = glmmTMB::nbinom1(), verbose = TRUE
# )
# summary(m2)
# res2 <- DHARMa::simulateResiduals(m2, plot = TRUE)

# continuous response distribution with dispersion varying by taxa? ------------
dpos <- d |>
  group_by(serial_number) |>
  mutate(any_zero = 0 %in% species_average_daily_view) |>
  filter(!any_zero) |>
  ungroup()
nrow(d) - nrow(dpos)

m1.3 <- glmmTMB(
  species_average_daily_view ~ 1 + celebrity +
    (1 + celebrity | taxonomic_group) + (1 | serial_number),
  dispformula = ~ 0 + taxonomic_group,
  data = dpos, family = Gamma(link = "log"), verbose = TRUE
)
m1.3$sdr
summary(m1.3)
b <- coef(m1.3)
b$cond$taxonomic_group
res1.3 <- DHARMa::simulateResiduals(m1.3, plot = TRUE)
b <- broom.mixed::tidy(m1.3, effects = "ran_vals")
filter(b, term == "celebrity") %>%
  mutate(estimate = estimate + fixef(m1.3)$cond[[2]]) %>%
  ggplot(aes(exp(estimate),
    y = level,
    xmin = exp(estimate - 2 * std.error), xmax = exp(estimate + 2 * std.error)
  )) +
  geom_pointrange() +
  geom_vline(xintercept = 1, lty = 2)

# what about for the upper 50% of celebrities? ------------------------

.med <- median(dpos$celeb_average_daily_view, na.rm = TRUE)
.mean <- mean(dpos$celeb_average_daily_view, na.rm = TRUE)
.quant80 <- as.numeric(quantile(dpos$celeb_average_daily_view, probs = 0.8, na.rm = TRUE))

dpos <- dpos |>
  group_by(serial_number) |>
  mutate(above_med = max(celeb_average_daily_view, na.rm = TRUE) >= .med ) |>
  mutate(above_mean = max(celeb_average_daily_view, na.rm = TRUE) >= .mean) |>
  mutate(above_80quant = max(celeb_average_daily_view, na.rm = TRUE) >= .quant80) |>
  ungroup()

dpos_med <- filter(dpos, above_med)
dpos_mean <- filter(dpos, above_mean)
dpos_quant <- filter(dpos, above_80quant)

nrow(dpos)
nrow(dpos_med)
nrow(dpos_mean)
nrow(dpos_quant)

m1.4 <- glmmTMB(
  species_average_daily_view ~ 1 + celebrity +
    (1 + celebrity | taxonomic_group) + (1 | serial_number),
  # dispformula = ~ 0 + taxonomic_group,
  data = dpos_med, family = Gamma(link = "log"), verbose = TRUE
)
m1.4$sdr
summary(m1.4)
b <- coef(m1.4)
b$cond$taxonomic_group
# res1.4 <- DHARMa::simulateResiduals(m1.4, plot = TRUE)
b <- broom.mixed::tidy(m1.4, effects = "ran_vals")
filter(b, term == "celebrity") %>%
  mutate(estimate = estimate + fixef(m1.4)$cond[[2]]) %>%
  ggplot(aes(exp(estimate),
    y = level,
    xmin = exp(estimate - 2 * std.error), xmax = exp(estimate + 2 * std.error)
  )) +
  geom_pointrange() +
  geom_vline(xintercept = 1, lty = 2)

# go brms? --------------------------------
library(brms)
# library(future)
# plan(multisession)
# fit1 <- brm(bf(symptom_post ~ group, sigma ~ group),
  # data = dat1, family = gaussian())

get_prior(bf(
  species_average_daily_view ~ 1 + celebrity +
    (1 + celebrity | taxonomic_group) + (1 | serial_number),
  shape ~ 0 + as.factor(taxonomic_group)),
  data = dpos, family = Gamma(link = "log"))

priors <-
  prior(normal(0, 1), class = "b") +
  prior(normal(0, 5), class = "Intercept") +
  prior(student_t(3, 0, 2.5), class = "sd") +
  prior(student_t(3, 0, 2.5), class = "b", dpar = "shape") +
  prior(lkj_corr_cholesky(1), class = "L")

fit_brms1.4 <- brm(
  bf(species_average_daily_view ~ 1 + celebrity +
      (1 + celebrity | taxonomic_group) + (1 | serial_number),
    shape ~ 0 + as.factor(taxonomic_group)),
  data = dpos,
  family = Gamma(link = "log"),
  backend = "cmdstanr",
  iter = 800L,
  chains = 6L,
  cores = future::availableCores(logical = FALSE),
  prior = priors
)
fit_brms1.4
# 3 divergences! look at shape parameter model/increase adapt_delta from 0.8

loo::loo(fit_brms1.3)
loo::loo(fit_brms1.4)

p1 <- tidybayes::spread_draws(fit_brms1.4, r_taxonomic_group[taxa, param])
p2 <- tidybayes::spread_draws(fit_brms1.4, b_celebrity)

p1_celeb <- filter(p1, param == "celebrity") |>
  left_join(p2)

ggplot(p1_celeb, aes(taxa, exp(b_celebrity + r_taxonomic_group))) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_violin(fill = NA) +
  coord_flip() + xlab("") + ylab("Rating") +
  theme_light()

group_by(p1_celeb, taxa) |>
  mutate(theta = exp(b_celebrity + r_taxonomic_group)) |>
  summarise(prob_gt_one = mean(theta > 1), CI95_lwr = quantile(theta, 0.025), median = quantile(theta, 0.5), CI95_upr = quantile(theta, 0.975)) |>
  knitr::kable(digits = 2L)

# continuous wiki counts: -----------------------------
d$log_celeb_views <- log(d$celeb_total_views)
d_celeb <- filter(d, celebrity == 1)

m3 <- glmmTMB(
  species_total_views ~ 1 + log_celeb_views +
    (1 + log_celeb_views | taxonomic_group),
  data = d_celeb, family = glmmTMB::nbinom2(), verbose = TRUE
)
summary(m3)
res <- DHARMa::simulateResiduals(m3, plot = TRUE)
b <- coef(m3)
b$cond$taxonomic_group

b <- broom.mixed::tidy(m3, effects = "ran_vals")
filter(b, term == "log_celeb_views") %>%
  mutate(estimate = estimate + fixef(m3)$cond[[2]]) %>%
  ggplot(aes(estimate,
    y = level,
    xmin = estimate - 2 * std.error, xmax = estimate + 2 * std.error
  )) +
  geom_pointrange() +
  geom_vline(xintercept = 0, lty = 2)

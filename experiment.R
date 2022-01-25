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

# paired; continuous wiki counts:
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

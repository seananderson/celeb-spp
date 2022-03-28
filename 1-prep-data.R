library(tidyverse)

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
names(d) <- tolower(names(d))
d$taxonomic_group <- as.factor(d$taxonomic_group)
dpos <- d |>
  group_by(serial_number) |>
  mutate(any_zero = 0 %in% species_average_daily_view) |>
  filter(!any_zero) |>
  ungroup()
nrow(d) - nrow(dpos)

.med <- median(dpos$celeb_average_daily_view, na.rm = TRUE)
.mean <- mean(dpos$celeb_average_daily_view, na.rm = TRUE)
.quant80 <- as.numeric(quantile(dpos$celeb_average_daily_view, probs = 0.8,
  na.rm = TRUE))

dpos <- dpos |>
  group_by(serial_number) |>
  mutate(above_med = max(celeb_average_daily_view, na.rm = TRUE) >= .med ) |>
  mutate(above_mean = max(celeb_average_daily_view, na.rm = TRUE) >= .mean) |>
  mutate(above_1000 = max(celeb_average_daily_view, na.rm = TRUE) >= 1000) |>
  mutate(above_80quant = max(celeb_average_daily_view, na.rm = TRUE) >= .quant80) |>
  ungroup()

dpos_med <- filter(dpos, above_med)
dpos_mean <- filter(dpos, above_mean)
dpos_quant <- filter(dpos, above_80quant)
dpos_1000 <- filter(dpos, above_1000)

nrow(dpos)
nrow(dpos_med)
nrow(dpos_mean)
nrow(dpos_1000)
nrow(dpos_quant)


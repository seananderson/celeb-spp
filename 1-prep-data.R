library(dplyr)

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

d <- group_by(d, serial_number) |>
  mutate(celeb_average_daily_view = mean(celeb_average_daily_view, na.rm = TRUE))

dpos <- d |>
  group_by(serial_number) |>
  mutate(any_zero = 0 %in% species_average_daily_view) |>
  filter(!any_zero) |>
  ungroup()
nrow(d) - nrow(dpos)

.med <- median(dpos$celeb_average_daily_view, na.rm = TRUE)
.mean <- mean(dpos$celeb_average_daily_view, na.rm = TRUE)
.quant80 <- as.numeric(quantile(dpos$celeb_average_daily_view,
  probs = 0.8,
  na.rm = TRUE
))

dpos <- dpos |>
  group_by(serial_number) |>
  mutate(above_med = max(celeb_average_daily_view, na.rm = TRUE) >= .med) |>
  mutate(above_mean = max(celeb_average_daily_view, na.rm = TRUE) >= .mean) |>
  mutate(above_1000 = max(celeb_average_daily_view, na.rm = TRUE) >= 1000) |>
  mutate(above_100 = max(celeb_average_daily_view, na.rm = TRUE) >= 100) |>
  mutate(above_10 = max(celeb_average_daily_view, na.rm = TRUE) >= 10) |>
  mutate(above_80quant = max(celeb_average_daily_view, na.rm = TRUE) >= .quant80) |>
  ungroup()

d <- d |>
  group_by(serial_number) |>
  mutate(above_1000 = max(celeb_average_daily_view, na.rm = TRUE) >= 1000) |>
  mutate(above_100 = max(celeb_average_daily_view, na.rm = TRUE) >= 100) |>
  mutate(above_10 = max(celeb_average_daily_view, na.rm = TRUE) >= 10) |>
  ungroup()

dpos_med <- filter(dpos, above_med)
dpos_mean <- filter(dpos, above_mean)
dpos_quant <- filter(dpos, above_80quant)
dpos_1000 <- filter(dpos, above_1000)
dpos_100 <- filter(dpos, above_100)
dpos_10 <- filter(dpos, above_10)

d_1000 <- filter(d, above_1000)
d_100 <- filter(d, above_100)
d_10 <- filter(d, above_10)

d_u1000 <- filter(d, !above_1000)
d_u100 <- filter(d, !above_100)
d_u10 <- filter(d, !above_10)

nrow(dpos)
nrow(dpos_med)
nrow(dpos_mean)
nrow(dpos_1000)
nrow(dpos_100)
nrow(dpos_10)
nrow(dpos_quant)


nrow(d_1000)
nrow(d_u100)
nrow(d_100)
nrow(d_10)
nrow(d_u10)

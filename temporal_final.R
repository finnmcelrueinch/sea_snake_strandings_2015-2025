#######################
## Temporal analysis ##
#######################
library(circular)
library(CircMLE)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ggplot2)
library(viridis)
library(scales)
library(ggrepel)
library(ggridges)

set.seed(123)


data_file <- "2026-02-23_cleaned_data.csv"

min_year <- 2015
max_year <- 2025

min_n_group_tables <- 20
min_n_group_plots  <- 10


max_gap_days <- 3
min_event_n  <- 3
mass_event_n <- 7
main_dist_km <- 50

master_data <- read_csv(data_file)

seasonal_data <- master_data %>%
  filter(
    date_source == "found",
    date_precision %in% c("day", "month"),
    !is.na(event_date),
    event_date >= as.Date(sprintf("%d-01-01", min_year)),
    event_date <= as.Date(sprintf("%d-12-31", max_year))
  ) %>%
  mutate(
    year = year(event_date),
    month_num = month(event_date),
    month = month(event_date, label = TRUE, abbr = FALSE),
    season = case_when(
      month_num %in% c(12, 1, 2)  ~ "Summer",
      month_num %in% c(3, 4, 5)   ~ "Autumn",
      month_num %in% c(6, 7, 8)   ~ "Winter",
      month_num %in% c(9, 10, 11) ~ "Spring"
    ),
    month = factor(month, levels = month.name)
  )


seasonal_data <- seasonal_data %>%
  mutate(
    event_date_approx = if_else(
      date_precision == "month",
      as.Date(sprintf("%04d-%02d-15", year(event_date), month(event_date))),
      event_date
    ),
    doy = yday(event_date_approx),
    days_in_year = if_else(leap_year(year(event_date_approx)), 366, 365),
    angle_rad = 2 * pi * (doy - 1) / days_in_year,
    angle = circular(angle_rad, units = "radians", template = "clock24")
  )


# Overall strandings
breaks_12 <- seq(0, 2 * pi, length.out = 13)

max_count <- max(
  hist(seasonal_data$angle_rad, breaks = breaks_12, plot = FALSE)$counts
)

p_rose_overall <- ggplot(seasonal_data, aes(x = angle_rad)) +
  geom_histogram(
    breaks = breaks_12,
    fill = "#1B9E77",
    color = "white",
    linewidth = 0.4
  ) +
  coord_polar(start = -pi / 12) +
  scale_x_continuous(
    breaks = breaks_12[-13],
    labels = month.abb
  ) +
  scale_y_continuous(
    breaks = pretty(c(0, max_count), n = 4),
    name = "Number of stranding records"
  ) +
  labs(
    title = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_rose_overall)

# Species-specific strandings

species_order <- species_month_props %>%
  group_by(taxon_label) %>%
  summarise(
    peak_month = month_num[which.max(prop)],
    peak_prop = max(prop),
    .groups = "drop"
  ) %>%
  arrange(peak_month, desc(peak_prop), taxon_label) %>%
  pull(taxon_label)

species_month_props <- species_month_props %>%
  mutate(taxon_label = factor(taxon_label, levels = rev(species_order)))


n_species <- n_distinct(species_month_props$taxon_label)
species_cols <- setNames(
  hcl.colors(n_species, palette = "Spectral"),
  levels(species_month_props$taxon_label)
)



species_month_smooth <- species_month_props %>%
  group_by(taxon_label) %>%
  group_modify(~{
    df <- .x %>% arrange(month_num)
    
    
    x_ext <- c(df$month_num - 12, df$month_num, df$month_num + 12)
    y_ext <- c(df$prop, df$prop, df$prop)
    
    
    spline_fit <- smooth.spline(x_ext, y_ext, spar = 0.45)
    
    x_new <- seq(1, 12, length.out = 240)
    y_new <- predict(spline_fit, x = x_new)$y
    y_new[y_new < 0] <- 0
    
    tibble(
      month_num = x_new,
      prop_smooth = y_new
    )
  }) %>%
  ungroup() %>%
  left_join(
    species_month_props %>%
      group_by(taxon_label) %>%
      summarise(
        peak_month = month_num[which.max(prop)],
        peak_prop = max(prop),
        .groups = "drop"
      ),
    by = "taxon_label"
  )


species_month_points <- species_month_props %>%
  group_by(taxon_label) %>%
  mutate(is_peak = prop == max(prop)) %>%
  ungroup()


p_species_prop <- ggplot(
  species_month_smooth,
  aes(
    x = month_num,
    y = taxon_label,
    height = prop_smooth,
    fill = taxon_label
  )
) +
  geom_ridgeline(
    stat = "identity",
    scale = 7,
    colour = "grey15",
    linewidth = 0.35,
    alpha = 0.95,
    min_height = 0
  ) +
  scale_fill_manual(values = species_cols, guide = "none") +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb,
    limits = c(1, 12),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Month",
    y = NULL,
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "transparent", linewidth = 0),
    axis.title.x = element_text(size = 13, colour = "black"),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black", face = "italic"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11),
    plot.margin = margin(10, 12, 10, 10)
  )

p_species_prop

ggsave(
  filename = "species_monthly_proportion_ridgeline.png",
  plot = p_species_prop,
  width = 8.27,
  height = 11.69,
  units = "in",
  dpi = 600,
  bg = "white"
)

#Hermans-Rasson test
wrap_angle_2pi <- function(angle_rad_numeric) {
  angle_rad_numeric %% (2 * pi)
}

angle_to_doy <- function(angle_rad_numeric) {
  angle_wrapped <- wrap_angle_2pi(angle_rad_numeric)
  (angle_wrapped / (2 * pi)) * 365.25 + 1
}

doy_to_date <- function(doy, origin_year = 2020) {
  as.Date(doy - 1, origin = sprintf("%d-01-01", origin_year))
}

angles_all <- circular(
  wrap_angle_2pi(as.numeric(seasonal_data$angle)),
  units = "radians",
  modulo = "2pi"
)

overall_r <- rho.circular(angles_all)
overall_hr <- HR_test(angles_all)
overall_hr_stat <- overall_hr[1]
overall_hr_p <- overall_hr[2]

cat("\n=== HR TEST FOR OVERALL SEASONALITY ===\n")
cat("Mean resultant length r =", overall_r, "\n")
cat("Hermans-Rasson test statistic =", overall_hr_stat, "\n")
cat("Hermans-Rasson p-value =", overall_hr_p, "\n")


circular_group_summary <- function(df, group_var, min_n = 1, origin_year = 2020) {
  g <- rlang::ensym(group_var)
  
  df %>%
    group_by(!!g) %>%
    reframe(
      n = n(),
      angle_grp = list(
        circular(
          wrap_angle_2pi(as.numeric(angle)),
          units = "radians",
          modulo = "2pi"
        )
      )
    ) %>%
    mutate(
      r = purrr::map_dbl(angle_grp, rho.circular),
      mean_angle_rad = purrr::map_dbl(angle_grp, ~ wrap_angle_2pi(as.numeric(mean.circular(.x)))),
      hr_stat = purrr::map_dbl(angle_grp, ~ HR_test(.x)[1]),
      hr_p = purrr::map_dbl(angle_grp, ~ HR_test(.x)[2])
    ) %>%
    dplyr::select(-angle_grp) %>%
    filter(n >= min_n) %>%
    mutate(
      hr_p_adj = p.adjust(hr_p, method = "BH"),
      peak_doy = angle_to_doy(mean_angle_rad),
      peak_date = doy_to_date(peak_doy, origin_year = origin_year)
    )
}

species_data <- seasonal_data %>%
  filter(taxon_label != "unidentified")

species_circ <- circular_group_summary(
  species_data,
  taxon_label,
  min_n = min_n_group_tables,
  origin_year = 2020
) %>%
  arrange(hr_p_adj, desc(r))

cat("\n=== SPECIES CIRCULAR SUMMARY (n >= ", min_n_group_tables, ") ===\n", sep = "")
print(species_circ, n = Inf)

# Interannual variation

monthly_by_year <- seasonal_data %>%
  count(year, month_num, name = "n") %>%
  tidyr::complete(
    year = full_seq(min_year:max_year, 1),
    month_num = 1:12,
    fill = list(n = 0)
  )

monthly_norm <- monthly_by_year %>%
  group_by(year) %>%
  mutate(
    total_year = sum(n),
    prop = if_else(total_year > 0, n / total_year, 0)
  ) %>%
  ungroup()



monthly_smooth <- monthly_norm %>%
  group_by(year) %>%
  group_modify(~{
    df <- .x %>% arrange(month_num)
    
    x_ext <- c(df$month_num - 12, df$month_num, df$month_num + 12)
    y_ext <- c(df$prop, df$prop, df$prop)
    
    spline_fit <- smooth.spline(x_ext, y_ext, spar = 0.45)
    
    x_new <- seq(1, 12, length.out = 240)
    y_new <- predict(spline_fit, x = x_new)$y
    y_new[y_new < 0] <- 0
    
    tibble(
      month_num = x_new,
      prop_smooth = y_new
    )
  }) %>%
  ungroup() %>%
  mutate(year_f = factor(year, levels = rev(sort(unique(year)))))

monthly_points <- monthly_norm %>%
  mutate(
    year_f = factor(year, levels = rev(sort(unique(year))))
  ) %>%
  group_by(year) %>%
  mutate(is_peak = prop == max(prop)) %>%
  ungroup()



p_ridge_year <- ggplot(
  monthly_smooth,
  aes(
    x = month_num,
    y = year_f,
    height = prop_smooth,
    fill = year_f
  )
) +
  geom_ridgeline(
    stat = "identity",
    scale = 8,
    colour = "grey15",
    linewidth = 0.35,
    alpha = 0.95,
    min_height = 0
  ) +
 scale_fill_viridis_d(option = "C", guide = "none") +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb,
    limits = c(1, 12),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Month",
    y = NULL
  ) +
  theme_ridges() +
  theme(
    text = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "transparent", linewidth = 0.3),
    axis.title.x = element_text(size = 13, colour = "black"),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11),
    plot.margin = margin(10, 12, 10, 10)
  )

print(p_ridge_year)


ggsave(
  filename = "interannual_ridgeplot.png",
  plot = p_ridge_year,
  width = 8.27,
  height = 11.69,
  units = "in",
  dpi = 600,
  bg = "white"
)



# Cluster Analysis
cluster_data <- seasonal_data %>%
  filter(
    !is.na(event_date_approx),
    !is.na(lat),
    !is.na(lon)
  ) %>%
  mutate(
    lat = as.numeric(lat),
    lon = as.numeric(lon)
  ) %>%
  arrange(event_date_approx) %>%
  mutate(record_id = row_number())

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  to_rad <- pi / 180
  
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  
  lat1 <- lat1 * to_rad
  lat2 <- lat2 * to_rad
  
  a <- sin(dlat / 2)^2 +
    cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  
  2 * R * asin(pmin(1, sqrt(a)))
}

n <- nrow(cluster_data)

date_num <- as.numeric(cluster_data$event_date_approx)
time_diff_mat <- abs(outer(date_num, date_num, "-"))

lat_mat1 <- matrix(cluster_data$lat, n, n)
lat_mat2 <- t(lat_mat1)
lon_mat1 <- matrix(cluster_data$lon, n, n)
lon_mat2 <- t(lon_mat1)

dist_mat <- haversine_km(lat_mat1, lon_mat1, lat_mat2, lon_mat2)

adj_mat <- (time_diff_mat <= max_gap_days) & (dist_mat <= main_dist_km)
diag(adj_mat) <- FALSE

component_id <- rep(NA_integer_, n)
current_component <- 0L

for (i in seq_len(n)) {
  if (is.na(component_id[i])) {
    current_component <- current_component + 1L
    queue <- i
    component_id[i] <- current_component
    
    while (length(queue) > 0) {
      v <- queue[1]
      queue <- queue[-1]
      neighbors <- which(adj_mat[v, ] & is.na(component_id))
      if (length(neighbors) > 0) {
        component_id[neighbors] <- current_component
        queue <- c(queue, neighbors)
      }
    }
  }
}

event_data_refined <- cluster_data %>%
  mutate(
    event_id = paste0("EV_", component_id),
    dist_threshold_km = main_dist_km
  )

event_summary_all <- event_data_refined %>%
  group_by(event_id) %>%
  summarise(
    start = min(event_date_approx),
    end = max(event_date_approx),
    duration = as.numeric(end - start) + 1,
    n_records = n(),
    centroid_lat = mean(lat, na.rm = TRUE),
    centroid_lon = mean(lon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    event_type = case_when(
      n_records >= mass_event_n & duration <= max_gap_days ~ "Mass stranding",
      n_records >= min_event_n ~ "Stranding event",
      TRUE ~ "Single stranding"
    )
  )

event_summary_refined <- event_summary_all %>%
  filter(n_records >= min_event_n)

mass_stranding_events <- event_summary_all %>%
  filter(n_records >= mass_event_n, duration <= max_gap_days)

cat("\n=== EVENT COUNTS ===\n")
cat("Stranding events (>= ", min_event_n, " records): ", nrow(event_summary_refined), "\n", sep = "")
cat("Mass stranding events (>= ", mass_event_n, " records within ", max_gap_days, " days): ",
    nrow(mass_stranding_events), "\n", sep = "")


event_summary_plot <- event_summary_all %>%
  mutate(
    event_type = case_when(
      n_records >= mass_event_n & duration <= max_gap_days ~ "Mass stranding",
      n_records >= min_event_n ~ "Stranding event",
      TRUE ~ "Single stranding"
    ),
    year_decimal = lubridate::year(start) +
      (lubridate::yday(start) - 1) / ifelse(lubridate::leap_year(start), 366, 365),
    date_lab = format(start, "%d %b %Y")
  )

single_plot <- event_summary_plot %>%
  filter(event_type == "Single stranding")

stranding_plot <- event_summary_plot %>%
  filter(event_type == "Stranding event")

mass_plot <- event_summary_plot %>%
  filter(event_type == "Mass stranding")

p_event_time <- ggplot() +
  geom_jitter(
    data = single_plot,
    aes(x = year_decimal, y = n_records, color = event_type),
    width = 0.08,
    height = 0.08,
    size = 2.0,
    alpha = 0.85
  ) +
  geom_jitter(
    data = stranding_plot,
    aes(x = year_decimal, y = n_records, color = event_type),
    width = 0.08,
    height = 0.08,
    size = 2.4,
    alpha = 0.9
  ) +
  geom_point(
    data = mass_plot,
    aes(x = year_decimal, y = n_records, color = event_type),
    shape = 21,
    fill = "#D95F02",
    size = 4.8,
    stroke = 0.8
  ) +
  geom_text_repel(
    data = mass_plot,
    aes(x = year_decimal, y = n_records, label = date_lab),
    size = 5,
    fontface = "bold",
    color = "#D95F02",
    box.padding = 0.55,
    point.padding = 0.4,
    force = 3,
    nudge_y = 0.35,
    max.overlaps = Inf,
    segment.color = NA,
    seed = 123
  ) +
  scale_color_manual(
    values = c(
      "Single stranding" = "grey75",
      "Stranding event" = "#2C7FB8",
      "Mass stranding" = "#D95F02"
    ),
    breaks = c("Single stranding", "Stranding event", "Mass stranding"),
    name = "Event type"
  ) +
  scale_x_continuous(
    breaks = seq(min_year, max_year, by = 1),
    labels = seq(min_year, max_year, by = 1),
    limits = c(min_year, max_year + 0.99),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq(
      floor(min(event_summary_plot$n_records)),
      ceiling(max(event_summary_plot$n_records)),
      by = 1
    ),
    labels = seq(
      floor(min(event_summary_plot$n_records)),
      ceiling(max(event_summary_plot$n_records)),
      by = 1
    )
  ) +
  labs(
    x = "Year",
    y = "Strandings per event",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_event_time)

ggsave(
  filename = "stranding_events.png",
  plot = p_event_time,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

# Index of Dispersion Test
daily_counts <- seasonal_data %>%
  count(event_date_approx, name = "n") %>%
  complete(
    event_date_approx = seq(min(event_date_approx), max(event_date_approx), by = "day"),
    fill = list(n = 0)
  ) %>%
  mutate(
    year = year(event_date_approx),
    doy = yday(event_date_approx)
  )

mean_daily <- mean(daily_counts$n)
var_daily <- var(daily_counts$n)
dispersion_index <- var_daily / mean_daily

n_days <- nrow(daily_counts)
chi_sq <- (n_days - 1) * dispersion_index
p_value <- pchisq(chi_sq, df = n_days - 1, lower.tail = FALSE)

cat("\n=== INDEX OF DISPERSION TEST ===\n")
cat("Index of dispersion =", round(dispersion_index, 3), "\n")
cat("Chi-square =", round(chi_sq, 3), "\n")

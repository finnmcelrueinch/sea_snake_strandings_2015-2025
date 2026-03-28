######################
## Spatial analysis ##
######################


library(tidyverse)
library(sf)
library(mapview)
library(terra)
library(spatstat)
library(ggspatial)
library(ozmaps)
library(prettymapr)
library(conflicted)
library(vegan)
library(raster)
library(ggplot2)
library(ggmap)
library(cowplot)
library(grid)
library(sfdep)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

clean_data <- read_csv("2026-02-23_cleaned_data.csv")

clean_sf <- 
  clean_data %>% 
  filter(!is.na(lon)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)


land <- ozmap_country


pop_rast <- rast("ABS_Population_Grid_2025.tif")
pop_rast <- project(pop_rast, "EPSG:4326")
pop_rast <- crop(pop_rast, vect(land))
pop_land <- mask(pop_rast, vect(land))
pop_df <- as.data.frame(pop_land, xy = TRUE, na.rm = TRUE)
names(pop_df)[3] <- "pop_density"

mapview(clean_sf, zcol = "taxon_label") + land

# Consistent plot theme
theme_thesis_map <- function() {
  theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 9, colour = "grey20"),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.25),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "#9fbfd4", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
      plot.margin = margin(t = 8, r = 2, b = 2, l = 2)
    )
}

sf_use_s2(F)

land_buffer <- 
  st_buffer(land, dist = 0.1) %>% 
  st_difference(land) %>% 
  st_cast("POLYGON") %>% 
  slice(-c(1,10, 11, 12, 17)) %>% 
  st_union()

hexbin <- 
  land_buffer %>% 
  st_make_grid(cellsize = 0.4, 
               flat_topped = T,
               square = F) %>% 
  st_as_sf() %>%
  mutate(a = st_intersects(., land_buffer) %>% as.numeric()) %>% 
  filter(a %in% 1) %>% 
  select(-a) %>% 
  st_transform(4326)

stranding_bins <-
  hexbin %>% 
  mutate(num_strandings = st_intersects(., clean_sf) %>% lengths()) %>% 
  rename(geometry = x)

state_boundaries <- st_boundary(ozmap_states)
coastline <- st_boundary(land)


# National hexbin map

bb <- st_bbox(land)

bb["ymin"] <- bb["ymin"] - 2
bb["ymax"] <- bb["ymax"] + 0.5
bb["xmin"] <- bb["xmin"] - 1
bb["xmax"] <- bb["xmax"] + 1

state_boundaries <- st_boundary(ozmap_states)

major_settlements <- tibble::tribble(
  ~city,        ~lon,      ~lat,
  "Darwin",     130.8456, -12.4634,
  "Cairns",     145.7781, -16.9186,
  "Townsville", 146.8179, -19.2589,
  "Brisbane",   153.0251, -27.4698,
  "Sydney",     151.2093, -33.8688,
  "Newcastle",  151.7817, -32.9283,
  "Melbourne",  144.9631, -37.8136,
  "Hobart",     147.3272, -42.8821,
  "Adelaide",   138.6007, -34.9285,
  "Perth",      115.8605, -31.9505
)

major_settlements_sf <- 
  st_as_sf(major_settlements, coords = c("lon", "lat"), crs = 4326) %>%
  mutate(feature = "Major Australian\nurban settlements")

theme_thesis_map <- function() {
  theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 12, colour = "grey20"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "#9fbfd4", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
    )
}

p_hexmap <- ggplot() +
  geom_sf(
    data = land,
    fill = "grey75",
    colour = "grey60",
    linewidth = 0.2
  ) +
  geom_sf(
    data = state_boundaries,
    colour = "grey65",
    linewidth = 0.18
  ) +
  layer_spatial(
    stranding_bins,
    aes(fill = num_strandings),
    col = "white",
    stroke = 0.12
  ) +
  geom_sf(
    data = major_settlements_sf,
    aes(shape = feature),
    size = 1,
    fill = "black",
    colour = "white",
    stroke = 0.3
  ) +
  scale_fill_distiller(
    palette = "Reds",
    direction = 1,
    trans = "log10",
    na.value = NA,
    name = "Total number of\nreported strandings",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(5, "cm"),
      barheight = unit(0.4, "cm")
    )
  ) +
  scale_shape_manual(
    values = c("Major Australian\nurban settlements" = 21),
    name = NULL
  ) +
  coord_sf(
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"]),
    expand = FALSE
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  annotation_scale(
    location = "br",
    width_hint = 0.16,
    text_cex = 1,
    line_width = 0.6
  ) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(0.9, "cm"),
    width = unit(0.9, "cm")
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        fill = "black",
        colour = "white",
        size = 4
      )
    )
  ) +
  theme_thesis_map()

p_hexmap


ggsave(
  filename = "stranding_hexmap.png",
  plot = p_hexmap,
  width = 8.27,
  height = 5.8,
  units = "in",
  dpi = 900,
  bg = "white"
)


# KDE analysis using MASS

library(MASS)

# Nationwide
st_kde <- function(points, cellsize, bandwidth, extent = NULL, expand = 0){
  require(MASS)
  require(raster)
  require(sf)
  
  if (is.null(extent)) {
    extent_vec <- st_bbox(points)[c(1, 3, 2, 4)]
  } else {
    extent_vec <- st_bbox(extent)[c(1, 3, 2, 4)]
  }
  

  extent_vec[1] <- extent_vec[1] - expand
  extent_vec[2] <- extent_vec[2] + expand
  extent_vec[3] <- extent_vec[3] - expand
  extent_vec[4] <- extent_vec[4] + expand
  
  n_y <- ceiling((extent_vec[4] - extent_vec[3]) / cellsize)
  n_x <- ceiling((extent_vec[2] - extent_vec[1]) / cellsize)
  
  extent_vec[2] <- extent_vec[1] + (n_x * cellsize) - cellsize
  extent_vec[4] <- extent_vec[3] + (n_y * cellsize) - cellsize
  
  coords <- st_coordinates(points)
  kde_mat <- kde2d(
    coords[,1],
    coords[,2],
    h = bandwidth,
    n = c(n_x, n_y),
    lims = extent_vec
  )
  
  raster(kde_mat)
}

points_kde <- st_kde(
  clean_sf,
  cellsize = 0.05,
  bandwidth = 1,
  extent = stranding_bins,
  expand = 6
)

kde_pol <- 
  points_kde %>% 
  rast() %>% 
  project("EPSG:4326") %>% 
  mask(
    vect(land_buffer %>% st_buffer(0.15)) %>% 
      project("EPSG:4326")
  )



kde_df <- terra::as.data.frame(kde_pol, xy = TRUE, na.rm = TRUE)
names(kde_df)[3] <- "kde"


bb <- st_bbox(land)
bb["xmin"] <- bb["xmin"] - 1
bb["xmax"] <- bb["xmax"] + 1
bb["ymin"] <- bb["ymin"] - 2
bb["ymax"] <- bb["ymax"] + 0.5


p_kde <- ggplot() +
  geom_raster(
    data = kde_df,
    aes(x = x, y = y, fill = kde, alpha = kde),
    interpolate = TRUE
  ) +
  scale_alpha_continuous(
    range = c(0.12, 1),
    guide = "none"
  ) +
  geom_sf(
    data = land,
    fill = "grey75",
    colour = "grey60",
    linewidth = 0.2
  ) +
  geom_sf(
    data = state_boundaries,
    colour = "grey65",
    linewidth = 0.18
  ) +
  scale_fill_gradientn(
    colours = c("navy", "blue", "cyan", "yellow", "orange", "red"),
    na.value = NA,
    name = "Kernel density value",
    breaks = c(0, 0.10, 0.20),
    labels = c("0", "0.10", "0.20"),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(5.5, "cm"),
      barheight = unit(0.45, "cm"),
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  scale_shape_manual(
    values = c("Major Australian\nurban settlements" = 21),
    name = NULL
  ) +
  coord_sf(
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"]),
    expand = FALSE
  ) +
  annotation_scale(
    location = "br",
    width_hint = 0.16,
    text_cex = 1,
    line_width = 0.6
  ) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(0.9, "cm"),
    width = unit(0.9, "cm")
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        fill = "black",
        colour = "white",
        size = 4
      )
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 12, colour = "grey20"),
    axis.ticks = element_line(colour = "grey30", linewidth = 0.25),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.spacing.y = unit(0.15, "cm"),
    legend.margin = margin(t = 2, b = 2),
    legend.box.margin = margin(b = 2),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "#a6c4d9", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
    plot.margin = margin(t = 8, r = 2, b = 2, l = 2)
  )

p_kde


ggsave(
  "nationwide_kde.png",
  p_kde,
  width = 8.27,
  height = 5.8,
  units = "in",
  dpi = 600,
  bg = "white"
)






#Per state


sf_use_s2(FALSE)
register_stadiamaps("956ecf85-d8a8-414b-a17a-ee98f8c68e17", write = FALSE)


panel_width  <- 4.13   
panel_height <- 5.85   
legend_width <- 4.13   
legend_height <- 1.25  
panel_dpi <- 600

clean_data <- 
  clean_data %>% 
  mutate(state = str_trim(str_to_upper(state)))

clean_sf <- 
  clean_data %>% 
  filter(!is.na(lon), !is.na(lat), !is.na(state)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

states_to_plot <- 
  clean_data %>% 
  filter(!is.na(state)) %>% 
  count(state, sort = TRUE) %>% 
  filter(n > 25) %>% 
  pull(state)

state_polys <- 
  ozmap_states %>% 
  st_transform(4326) %>% 
  mutate(
    state = recode(
      NAME,
      "Queensland" = "QLD",
      "New South Wales" = "NSW",
      "Victoria" = "VIC",
      "Tasmania" = "TAS",
      "South Australia" = "SA",
      "Western Australia" = "WA",
      "Northern Territory" = "NT",
      "Australian Capital Territory" = "ACT"
    )
  ) %>% 
  filter(state %in% states_to_plot) %>% 
  st_make_valid()

points_by_state <- 
  clean_sf %>% 
  filter(state %in% states_to_plot) %>% 
  split(.$state)



pop_rast <- rast("ABS_Population_Grid_2025.tif")
pop_rast <- project(pop_rast, "EPSG:4326")


st_kde <- function(points, cellsize, bandwith, extent = NULL) {
  require(MASS)
  require(raster)
  require(sf)
  
  if (is.null(extent)) {
    extent_vec <- st_bbox(points)[c(1, 3, 2, 4)]
  } else {
    extent_vec <- st_bbox(extent)[c(1, 3, 2, 4)]
  }
  
  n_y <- ceiling((extent_vec[4] - extent_vec[3]) / cellsize)
  n_x <- ceiling((extent_vec[2] - extent_vec[1]) / cellsize)
  
  extent_vec[2] <- extent_vec[1] + (n_x * cellsize) - cellsize
  extent_vec[4] <- extent_vec[3] + (n_y * cellsize) - cellsize
  
  coords <- st_coordinates(points)
  kde_mat <- kde2d(
    coords[, 1],
    coords[, 2],
    h = bandwith,
    n = c(n_x, n_y),
    lims = extent_vec
  )
  
  raster(kde_mat)
}

get_kde_params <- function(n) {
  case_when(
    n >= 200 ~ list(cellsize = 0.03, bandwith = 0.6),
    n >= 100 ~ list(cellsize = 0.04, bandwith = 0.8),
    TRUE     ~ list(cellsize = 0.05, bandwith = 1.0)
  )
}


theme_panel_map <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 12, colour = "grey20"),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.25),
      legend.position = "none",
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = NA, colour = NA),
      plot.margin = margin(4, 4, 4, 4)
    )
}

theme_kde_panel <- function() {
  theme_panel_map() +
    theme(
      panel.background = element_rect(fill = "#a6c4d9", colour = NA)
    )
}

theme_legend_only <- function() {
  theme_void(base_size = 11) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
}



make_raw_legend_plot <- function() {
  ggplot(
    tibble(x = 1, y = 1, grp = "Stranding record"),
    aes(x = x, y = y, colour = grp)
  ) +
    geom_point(size = 2.5, alpha = 1) +
    scale_colour_manual(
      values = c("Stranding record" = "black"),
      name = NULL
    ) +
    theme_legend_only()
}

make_population_legend_plot <- function() {
  pop_break_vals <- c(1, 10, 100, 1000, 10000)
  pop_breaks <- log10(pop_break_vals)
  
  ggplot(
    data.frame(
      x = seq_along(pop_breaks),
      y = 1,
      z = seq(min(pop_breaks), max(pop_breaks), length.out = length(pop_breaks))
    )
  ) +
    geom_tile(aes(x = x, y = y, fill = z)) +
    scale_fill_gradientn(
      colours = c("black", "purple", "deeppink", "orange", "yellow", "white"),
      limits = range(pop_breaks),
      breaks = pop_breaks,
      labels = scales::label_number(big.mark = ",")(pop_break_vals),
      name = "Population density (people/km²)",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(7, "cm"),
        barheight = unit(0.6, "cm"),
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    theme_legend_only()
}

make_kde_legend_plot <- function(kde_limits) {
  kde_breaks <- seq(kde_limits[1], kde_limits[2], length.out = 5)
  
  ggplot(
    data.frame(
      x = seq_along(kde_breaks),
      y = 1,
      z = seq(kde_limits[1], kde_limits[2], length.out = length(kde_breaks))
    )
  ) +
    geom_tile(aes(x = x, y = y, fill = z)) +
    scale_fill_gradientn(
      colours = c("navy", "blue", "cyan", "yellow", "orange", "red"),
      limits = kde_limits,
      breaks = kde_breaks,
      labels = scales::label_number(accuracy = 0.01)(kde_breaks),
      name = "Kernel density value",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(7, "cm"),
        barheight = unit(0.6, "cm"),
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    theme_legend_only()
}



make_state_panels <- function(pts_state, st, pop_rast) {

  
  poly_state <- 
    state_polys %>% 
    filter(state == st) %>% 
    st_transform(3577) %>% 
    st_union() %>% 
    st_make_valid() %>% 
    st_transform(4326)
  
  n_pts <- nrow(pts_state)
  pars <- get_kde_params(n_pts)
  
  
  poly_state_prj <- st_transform(poly_state, 3577)
  land_prj <- st_transform(ozmap_country, 3577)
  
 
  
  plot_extent_prj <- st_buffer(poly_state_prj, dist = 100000)
  plot_extent <- st_transform(plot_extent_prj, 4326)
  bb <- st_bbox(plot_extent)
  

  
  bbox_vec <- c(
    left   = unname(bb["xmin"]),
    bottom = unname(bb["ymin"]),
    right  = unname(bb["xmax"]),
    top    = unname(bb["ymax"])
  )
  
  raw_basemap <- get_stadiamap(
    bbox = bbox_vec,
    zoom = 6,
    maptype = "stamen_terrain",
    crop = TRUE,
    messaging = FALSE
  )
  
  pop_basemap <- get_stadiamap(
    bbox = bbox_vec,
    zoom = 6,
    maptype = "alidade_smooth_dark",
    crop = TRUE,
    messaging = FALSE
  )
  
  
  
  x_breaks <- pretty(c(bb["xmin"], bb["xmax"]), n = 4)
  y_breaks <- pretty(c(bb["ymin"], bb["ymax"]), n = 5)
  
  lon_label <- function(x) paste0(round(x, 0), "°E")
  lat_label <- function(y) paste0(abs(round(y, 0)), "°S")
  
  
  
  x_range <- as.numeric(bb["xmax"] - bb["xmin"])
  y_range <- as.numeric(bb["ymax"] - bb["ymin"])
  
  if (st %in% c("QLD", "NSW")) {
    sb_xmin <- bb["xmin"] + 0.04 * x_range
    sb_xmax <- bb["xmin"] + 0.28 * x_range
  } else {
    sb_xmin <- bb["xmax"] - 0.28 * x_range
    sb_xmax <- bb["xmax"] - 0.04 * x_range
  }
  
  sb_ymin <- bb["ymin"] + 0.03 * y_range
  sb_ymax <- bb["ymin"] + 0.05 * y_range
  sb_breaks <- seq(sb_xmin, sb_xmax, length.out = 5)
  
  scale_bar_layers <- list(
    annotate(
      "rect",
      xmin = sb_xmin - 0.01 * x_range,
      xmax = sb_xmax + 0.01 * x_range,
      ymin = sb_ymin - 0.015 * y_range,
      ymax = sb_ymax + 0.035 * y_range,
      fill = alpha("white", 0.75),
      colour = NA,
      linewidth = 0.2
    ),
    annotate("rect", xmin = sb_breaks[1], xmax = sb_breaks[2], ymin = sb_ymin, ymax = sb_ymax,
             fill = "black", colour = "black"),
    annotate("rect", xmin = sb_breaks[2], xmax = sb_breaks[3], ymin = sb_ymin, ymax = sb_ymax,
             fill = "white", colour = "black"),
    annotate("rect", xmin = sb_breaks[3], xmax = sb_breaks[4], ymin = sb_ymin, ymax = sb_ymax,
             fill = "black", colour = "black"),
    annotate("rect", xmin = sb_breaks[4], xmax = sb_breaks[5], ymin = sb_ymin, ymax = sb_ymax,
             fill = "white", colour = "black"),
    annotate(
      "text",
      x = sb_xmin,
      y = sb_ymax + 0.015 * y_range,
      label = "300 km",
      hjust = 0,
      size = 3,
      colour = "black"
    )
  )
  
  
  
  raw_df <- pts_state %>% 
    st_coordinates() %>% 
    as_tibble() %>% 
    rename(x = X, y = Y)
  
 
  
  pop_state <- crop(pop_rast, vect(plot_extent))
  pop_state <- mask(pop_state, vect(poly_state))
  
  pop_df <- terra::as.data.frame(pop_state, xy = TRUE, na.rm = TRUE)
  names(pop_df)[3] <- "population"
  
  pop_df <- pop_df %>% 
    mutate(
      population = as.numeric(population),
      pop_plot = if_else(population > 0, log10(population), NA_real_)
    )
  
  pop_res <- terra::res(pop_state)
  pop_w <- pop_res[1]
  pop_h <- pop_res[2]
  
  pop_break_vals <- c(1, 10, 100, 1000, 10000)
  pop_breaks <- log10(pop_break_vals)
  pop_limits <- range(pop_breaks)
  

  
  kde_extent <- st_as_sfc(st_bbox(plot_extent))
  
  kde_rast <- 
    st_kde(
      points = pts_state,
      cellsize = pars$cellsize,
      bandwith = pars$bandwith,
      extent = kde_extent
    ) %>% 
    rast()
  
  crs(kde_rast) <- "EPSG:4326"
  
  marine_band_prj <- 
    st_buffer(land_prj, dist = 90000) %>% 
    st_difference(land_prj) %>% 
    st_make_valid()
  
  coastal_mask_prj <- 
    st_intersection(
      marine_band_prj,
      st_buffer(poly_state_prj, dist = 60000)
    ) %>% 
    st_make_valid()
  
  coastal_mask <- st_transform(coastal_mask_prj, 4326)
  
  kde_offshore <- mask(kde_rast, vect(coastal_mask))
  
  kde_df <- terra::as.data.frame(kde_offshore, xy = TRUE, na.rm = TRUE)
  names(kde_df)[3] <- "kde"
  kde_df$kde <- as.numeric(kde_df$kde)
  
  kde_res <- terra::res(kde_offshore)
  kde_w <- kde_res[1]
  kde_h <- kde_res[2]
  
  kde_limits <- range(kde_df$kde, na.rm = TRUE)
  kde_breaks <- seq(kde_limits[1], kde_limits[2], length.out = 5)
  

  
  p_raw <- ggplot() +
    annotation_raster(
      raw_basemap,
      xmin = bb["xmin"], xmax = bb["xmax"],
      ymin = bb["ymin"], ymax = bb["ymax"]
    ) +
    geom_point(
      data = raw_df,
      aes(x = x, y = y),
      size = 1.0,
      alpha = 0.9,
      colour = "black"
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = lon_label,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = lat_label,
      expand = c(0, 0)
    ) +
    scale_bar_layers +
    coord_quickmap(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    ) +
    theme_panel_map()
  

  
  p_pop <- ggplot() +
    annotation_raster(
      pop_basemap,
      xmin = bb["xmin"], xmax = bb["xmax"],
      ymin = bb["ymin"], ymax = bb["ymax"]
    ) +
    geom_tile(
      data = pop_df,
      aes(x = x, y = y, fill = pop_plot),
      width = pop_w,
      height = pop_h,
      alpha = 0.75
    ) +
    scale_fill_gradientn(
      colours = c("black", "purple", "deeppink", "orange", "yellow", "white"),
      limits = pop_limits,
      breaks = pop_breaks,
      labels = scales::label_number(big.mark = ",")(pop_break_vals),
      na.value = NA,
      name = "Population density/km"
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = lon_label,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = lat_label,
      expand = c(0, 0)
    ) +
    coord_quickmap(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    ) +
    theme_panel_map()
  

  p_kde <- ggplot() +
    geom_tile(
      data = kde_df,
      aes(x = x, y = y, fill = kde),
      width = kde_w,
      height = kde_h,
      alpha = 0.9
    ) +
    geom_sf(
      data = poly_state,
      fill = "grey75",
      colour = "grey55",
      linewidth = 0.3
    ) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    ) +
    scale_fill_gradientn(
      name= "Kernel Density Value",
      colours = c("navy", "blue", "cyan", "yellow", "orange", "red"),
      limits = kde_limits,
      breaks = kde_breaks,
      labels = scales::label_number(accuracy = 0.01)(kde_breaks),
      na.value = NA
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = lon_label,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = lat_label,
      expand = c(0, 0)
    ) +
    theme_kde_panel()

  
  leg_raw <- make_raw_legend_plot()
  leg_pop <- make_population_legend_plot()
  leg_kde <- make_kde_legend_plot(kde_limits)
  
  list(
    raw = p_raw,
    pop = p_pop,
    kde = p_kde,
    legend_raw = leg_raw,
    legend_pop = leg_pop,
    legend_kde = leg_kde
  )
}



state_panels <- imap(points_by_state, ~ make_state_panels(.x, .y, pop_rast))

iwalk(state_panels, function(p_list, st) {
  
  ggsave(
    filename = paste0("state_raw_", st, ".png"),
    plot = p_list$raw,
    width = 5,
    height = 5.85,
    units = "in",
    dpi = panel_dpi,
    bg = "white"
  )
  
  ggsave(
    filename = paste0("state_population_", st, ".png"),
    plot = p_list$pop,
    width = 5,
    height = 5.85,
    units = "in",
    dpi = panel_dpi,
    bg = "white"
  )
  
  ggsave(
    filename = paste0("state_kde_", st, ".png"),
    plot = p_list$kde,
    width = 5,
    height = 5.85,
    units = "in",
    dpi = panel_dpi,
    bg = "white"
  )
})


iwalk(state_panels, function(p_list, st) {
  
  ggsave(
    filename = paste0("legend_raw_", st, ".png"),
    plot = p_list$legend_raw,
    width = legend_width,
    height = legend_height,
    units = "in",
    dpi = panel_dpi,
    bg = "white"
  )
  
  ggsave(
    filename = paste0("legend_population_", st, ".png"),
    plot = p_list$legend_pop,
    width = legend_width,
    height = legend_height,
    units = "in",
    dpi = panel_dpi,
    bg = "white"
  )
  
  ggsave(
    filename = paste0("legend_kde_", st, ".png"),
    plot = p_list$legend_kde,
    width = legend_width,
    height = legend_height,
    units = "in",
    dpi = panel_dpi,
    bg = "white"
  )
})


# Hot spot analysis (calculating the Gi* metric)

mapview(stranding_bins)

stranding_nb <- 
  stranding_bins %>% 
  mutate(nb = include_self(st_contiguity(geometry)),
         wt = st_weights(nb),
         stranding_lag = st_lag(num_strandings, nb, wt))

stranding_nb %>% 
  ggplot() +
  annotation_map_tile('cartolight', zoom = 5) +
  geom_sf(aes(fill = stranding_lag))

global_g_test(x = stranding_nb$num_strandings, 
              nb = stranding_nb$nb, 
              wt = stranding_nb$wt)


stranding_hot_spots <-
  stranding_nb %>% 
  mutate(
    gistar = local_g_perm(num_strandings, nb, wt, nsim = 499)
  ) %>% 
  unnest(gistar)

stranding_hot_spots %>% 
  ggplot() +
  annotation_map_tile('cartolight', zoom = 5) +
  geom_sf(aes(fill = gi)) +
  scale_fill_gradient2()

# very (hot/cold): p < 0.01
# hot/cold: p <= 0.05
# somewhat (hot/cold): p <= 0.1

strandings_classified <-
  stranding_hot_spots %>% 
  dplyr::select(gi, p_folded_sim) %>% 
  mutate(
    classification = 
      case_when(
        gi > 0 & p_folded_sim <= 0.01 ~ "Very hot (99%)",
        gi > 0 & p_folded_sim <= 0.05 ~ "Hot (95%)",
        gi > 0 & p_folded_sim <= 0.1  ~ "Somewhat hot (90%)",
        gi < 0 & p_folded_sim <= 0.01 ~ "Very cold (99%)",
        gi < 0 & p_folded_sim <= 0.05 ~ "Cold (95%)",
        gi < 0 & p_folded_sim <= 0.1  ~ "Somewhat cold (90%)",
        TRUE ~ "Insignificant"
      ),
    classification =
      factor(
        classification,
        levels = c(
          "Very hot (99%)",
          "Hot (95%)",
          "Somewhat hot (90%)",
          "Insignificant",
          "Somewhat cold (90%)",
          "Cold (95%)",
          "Very cold (99%)"
        )
      )
  )


# Plot classified Gi* hotspots - national map, styled consistently with other maps

state_boundaries <- 
  ozmap_states %>% 
  st_transform(4326) %>% 
  st_boundary()

bb <- st_bbox(land)
bb["xmin"] <- bb["xmin"] - 1
bb["xmax"] <- bb["xmax"] + 1
bb["ymin"] <- bb["ymin"] - 2
bb["ymax"] <- bb["ymax"] + 0.5

p_gistar_national <- ggplot() +
  geom_sf(
    data = land,
    fill = "grey75",
    colour = "grey60",
    linewidth = 0.2
  ) +
  geom_sf(
    data = state_boundaries,
    colour = "grey65",
    linewidth = 0.18
  ) +
  geom_sf(
    data = strandings_classified,
    aes(fill = classification),
    colour = "white",
    linewidth = 0.08,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "Very hot (99%)" = "#b2182b",
      "Hot (95%)" = "#ef8a62",
      "Somewhat hot (90%)" = "#fddbc7",
      "Insignificant" = "grey90",
      "Somewhat cold (90%)" = "#d1e5f0",
      "Cold (95%)" = "#67a9cf",
      "Very cold (99%)" = "#2166ac"
    ),
    drop = FALSE,
    name = "Gi* hotspot classification"
  ) +
  coord_sf(
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"]),
    expand = FALSE
  ) +
  annotation_scale(
    location = "br",
    width_hint = 0.16,
    text_cex = 1,
    line_width = 0.6
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(
    title = NULL,
    subtitle = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 12, colour = "grey20"),
    axis.ticks = element_line(colour = "grey30", linewidth = 0.25),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "#9fbfd4", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
    plot.margin = margin(t = 8, r = 2, b = 2, l = 2)
  )

p_gistar_national

ggsave(
  "nationwide_Gi_star.png",
  p_gistar_national,
  width =8,
  height = 5.8,
  units = "in",
  dpi = 600,
  bg = "white"
)


# Per-state Gi* hotspot maps

sf_use_s2(FALSE)


clean_data <- 
  clean_data %>% 
  mutate(state = str_trim(str_to_upper(state)))

clean_sf <- 
  clean_data %>% 
  filter(!is.na(lon), !is.na(lat), !is.na(state)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)


states_gi <- 
  clean_data %>% 
  filter(!is.na(state)) %>% 
  count(state, sort = TRUE) %>% 
  filter(n >= 25) %>% 
  pull(state)



state_polys <- 
  ozmap_states %>% 
  st_transform(4326) %>% 
  mutate(
    state = recode(
      NAME,
      "Queensland" = "QLD",
      "New South Wales" = "NSW",
      "Victoria" = "VIC",
      "Tasmania" = "TAS",
      "South Australia" = "SA",
      "Western Australia" = "WA",
      "Northern Territory" = "NT",
      "Australian Capital Territory" = "ACT"
    )
  ) %>% 
  filter(state %in% states_gi) %>% 
  st_make_valid()



points_by_state_gi <- 
  clean_sf %>% 
  filter(state %in% states_gi) %>% 
  split(.$state)



get_hex_params <- function(n) {
  dplyr::case_when(
    n >= 200 ~ list(cellsize = 15000),  
    n >= 100 ~ list(cellsize = 20000),  
    TRUE     ~ list(cellsize = 25000)   
  )
}



lon_label <- function(x) paste0(round(x, 0), "°E")
lat_label <- function(y) paste0(abs(round(y, 0)), "°S")



theme_gi_state <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 10, colour = "grey20"),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.25),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "#9fbfd4", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
      plot.margin = margin(t = 8, r = 2, b = 2, l = 2)
    )
}



gi_levels <- c(
  "Very hot (99%)",
  "Hot (95%)",
  "Somewhat hot (90%)",
  "Insignificant",
  "Somewhat cold (90%)",
  "Cold (95%)",
  "Very cold (99%)"
)

gi_cols <- c(
  "Very hot (99%)" = "#b2182b",
  "Hot (95%)" = "#ef8a62",
  "Somewhat hot (90%)" = "#fddbc7",
  "Insignificant" = "grey90",
  "Somewhat cold (90%)" = "#d1e5f0",
  "Cold (95%)" = "#67a9cf",
  "Very cold (99%)" = "#2166ac"
)


get_state_axis_breaks <- function(st, bb) {
  
  if (st == "QLD") {
    x_all <- seq(140, 150, by = 5)
    y_all <- seq(-30, -10, by = 5)
    
  } else if (st == "NSW") {
    x_all <- c(140, 145, 150)
    y_all <- seq(-38, -28, by = 2)
    
  } else if (st == "WA") {
    x_all <- seq(115, 130, by = 5)
    y_all <- seq(-35, -15, by = 5)
    
  } else if (st == "NT") {
    x_all <- seq(128, 138, by = 2)
    y_all <- seq(-25, -15, by = 5)
    
  } else {
    x_all <- pretty(c(bb["xmin"], bb["xmax"]), n = 4)
    y_all <- pretty(c(bb["ymin"], bb["ymax"]), n = 5)
  }
  
  x_breaks <- x_all[x_all >= bb["xmin"] & x_all <= bb["xmax"]]
  y_breaks <- y_all[y_all >= bb["ymin"] & y_all <= bb["ymax"]]
  
  list(x = x_breaks, y = y_breaks)
}


state_gistar_plots <- purrr::imap(points_by_state_gi, function(pts_state, st) {
  
  n_pts <- nrow(pts_state)
  pars <- get_hex_params(n_pts)
  

  pts_state_prj <- 
    pts_state %>% 
    st_transform(3577)
  

  poly_state_prj <- 
    state_polys %>% 
    filter(state == st) %>% 
    st_transform(3577) %>% 
    st_union() %>% 
    st_make_valid() %>% 
    st_collection_extract("POLYGON")
  

  land_prj <- 
    ozmap_country %>% 
    st_transform(3577) %>% 
    st_make_valid()
  
  marine_band_prj <- 
    st_buffer(land_prj, dist = 20000) %>% 
    st_difference(land_prj) %>% 
    st_make_valid()
  
  offshore_zone <- 
    st_intersection(marine_band_prj, poly_state_prj) %>% 
    st_make_valid() %>% 
    filter(!st_is_empty(.))
  

  hex_grid <- st_make_grid(
    offshore_zone,
    cellsize = pars$cellsize,
    flat_topped = TRUE,
    square = FALSE
  )
  
  hexbin_state <- 
    st_sf(geometry = hex_grid) %>% 
    mutate(a = lengths(st_intersects(geometry, offshore_zone))) %>% 
    filter(a > 0) %>% 
    dplyr::select(-a) %>% 
    st_make_valid() %>% 
    filter(!st_is_empty(.))
  

  stranding_bins_state <- 
    hexbin_state %>% 
    mutate(num_strandings = lengths(st_intersects(geometry, pts_state_prj)))
  

  stranding_nb_state <- 
    stranding_bins_state %>% 
    mutate(
      nb = include_self(st_contiguity(geometry)),
      wt = st_weights(nb),
      stranding_lag = st_lag(num_strandings, nb, wt)
    )
  

  stranding_hot_spots_state <- 
    stranding_nb_state %>% 
    mutate(
      gistar = local_g_perm(num_strandings, nb, wt, nsim = 499)
    ) %>% 
    tidyr::unnest(gistar)
  

  strandings_classified_plot <- 
    stranding_hot_spots_state %>% 
    dplyr::select(geometry, gi, p_folded_sim) %>% 
    mutate(
      classification = case_when(
        gi > 0 & p_folded_sim <= 0.01 ~ "Very hot (99%)",
        gi > 0 & p_folded_sim <= 0.05 ~ "Hot (95%)",
        gi > 0 & p_folded_sim <= 0.10 ~ "Somewhat hot (90%)",
        gi < 0 & p_folded_sim <= 0.01 ~ "Very cold (99%)",
        gi < 0 & p_folded_sim <= 0.05 ~ "Cold (95%)",
        gi < 0 & p_folded_sim <= 0.10 ~ "Somewhat cold (90%)",
        TRUE ~ "Insignificant"
      ),
      classification = factor(classification, levels = gi_levels)
    ) %>% 
    st_transform(4326)
  

  poly_state_plot <- 
    state_polys %>% 
    filter(state == st)
  

  poly_state_extent_prj <- st_transform(poly_state_plot, 3577)
  hotspot_extent_prj <- st_transform(strandings_classified_plot, 3577)
  
  plot_extent_prj <- c(
    st_geometry(poly_state_extent_prj),
    st_geometry(hotspot_extent_prj)
  ) %>% 
    st_union() %>% 
    st_buffer(40000)
  
  plot_extent <- st_transform(plot_extent_prj, 4326)
  bb <- st_bbox(plot_extent)
  

  ax <- get_state_axis_breaks(st, bb)
  

  ggplot() +
    geom_sf(
      data = poly_state_plot,
      fill = "grey75",
      colour = "grey60",
      linewidth = 0.25
    ) +
    geom_sf(
      data = strandings_classified_plot,
      aes(fill = classification),
      colour = "white",
      linewidth = 0.04,
      show.legend = TRUE
    ) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE,
      clip = "on"
    ) +
    scale_x_continuous(
      breaks = ax$x,
      labels = lon_label,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = ax$y,
      labels = lat_label,
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = gi_cols,
      drop = FALSE,
      name = "Gi* hotspot classification"
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(
      title = NULL,
      subtitle = NULL
    ) +
    theme_gi_state()
})



purrr::walk(state_gistar_plots, print)



panel_width <- 8
panel_height <- 5
panel_dpi <- 600

purrr::iwalk(state_gistar_plots, function(p, st) {
  ggsave(
    filename = paste0("state_gi_", st, ".png"),
    plot = p,
    width = panel_width,
    height = panel_height,
    units = "in",
    dpi = panel_dpi,
    bg = "white"
  )
})



legend_gi_plot <- ggplot(
  tibble(
    classification = factor(gi_levels, levels = gi_levels),
    x = 1,
    y = seq_along(gi_levels)
  ),
  aes(x = x, y = y, fill = classification)
) +
  geom_tile() +
  scale_fill_manual(
    values = gi_cols,
    drop = FALSE,
    name = "Gi* hotspot\nclassification"
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(colour = NA)
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.box.margin = margin(2, 2, 2, 2)
  )

legend_gi <- cowplot::get_legend(legend_gi_plot)

grid::grid.newpage()
grid::grid.draw(legend_gi)

png(
  filename = "legend_gi_vertical.png",
  width = 1200,
  height = 2200,
  res = 600,
  bg = "transparent"
)
grid::grid.newpage()
grid::grid.draw(legend_gi)
dev.off()


# Community composition of strandings by state

library(tidyverse)
library(vegan)
library(ggplot2)


comm_data <- 
  clean_data %>% 
  mutate(
    state = str_trim(str_to_upper(state)),
    taxon_label = str_trim(taxon_label)
  ) %>% 
  filter(
    !is.na(state),
    !is.na(taxon_label),
    
    
    !str_detect(str_to_lower(taxon_label), "unidentified"),
    
    
    !str_detect(taxon_label, "\\ssp\\.?$")
  )



states_keep <- 
  clean_data %>% 
  filter(!is.na(state)) %>% 
  count(state) %>% 
  filter(n >= 25) %>% 
  pull(state)

states_keep

state_species_filtered <- 
  clean_data %>%
  filter(state %in% states_keep) %>%
  filter(!str_detect(taxon_label, "sp\\.")) %>%   
  filter(!str_detect(taxon_label, "unidentified")) %>%
  count(state, taxon_label) %>%
  pivot_wider(
    names_from = taxon_label,
    values_from = n,
    values_fill = 0
  )


sort(unique(comm_data$taxon_label))

comm_matrix_filtered <- 
  state_species_filtered %>%
  column_to_rownames("state") %>%
  as.matrix()

bray_dist_filtered <- vegdist(comm_matrix_filtered, method = "bray")

as.matrix(bray_dist_filtered)



bray_matrix_filtered <- as.matrix(bray_dist_filtered)

bray_long_filtered <- 
  as.data.frame(as.table(bray_matrix_filtered)) %>%
  rename(
    state1 = Var1,
    state2 = Var2,
    dissimilarity = Freq
  )

ggplot(bray_long_filtered, aes(x = state1, y = state2, fill = dissimilarity)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(
    option = "C",
    name = "Bray-Curtis\nDissimilarity",
    limits = c(0, 1)
  ) +
  labs(
    x = "State",
    y = "State",
    title = "Bray-Curtis dissimilarity of stranded sea snake assemblages",
    subtitle = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  where <- match(flag, args)
  if (is.na(where) || where == length(args)) {
    stop(paste("missing", flag), call. = FALSE)
  }
  args[[where + 1L]]
}
optional_value_after <- function(flag) {
  where <- match(flag, args)
  if (is.na(where)) return(NA_character_)
  if (where == length(args)) stop(paste("missing", flag), call. = FALSE)
  args[[where + 1L]]
}

run_directory <- normalizePath(value_after("--run-directory"), mustWork = TRUE)
output <- value_after("--output")
configured_path <- optional_value_after("--config")
zoom_text <- optional_value_after("--zoom-radius")
zoom_radius <- if (is.na(zoom_text)) NA_real_ else as.numeric(zoom_text)
if (!is.na(zoom_radius) && (!is.finite(zoom_radius) || zoom_radius <= 0.0)) {
  stop("--zoom-radius must be a positive number", call. = FALSE)
}
field_paths <- sort(list.files(
  file.path(run_directory, "fields"),
  pattern = "^field_[0-9]+[.]csv$",
  full.names = TRUE
))
if (length(field_paths) == 0L) stop("no field snapshots found", call. = FALSE)
field <- fread(field_paths[[length(field_paths)]])
config_path <- if (is.na(configured_path)) {
  file.path(run_directory, "config", "continuum_requested.yaml")
} else {
  configured_path
}
config <- read_yaml(config_path)
activated_multiplier <- as.numeric(
  config$migration$activated_r_mobility_multiplier
)

nx <- as.integer(config$grid$shape[[1]])
ny <- as.integer(config$grid$shape[[2]])
nz <- as.integer(config$grid$shape[[3]])
if (nz != 1L) stop("this renderer expects a one-plane field", call. = FALSE)
origin_x <- as.numeric(config$grid$origin[[1]])
origin_y <- as.numeric(config$grid$origin[[2]])
spacing <- as.numeric(config$grid$spacing_voxels)
stride <- as.integer(config$output$field_stride)
x_values <- origin_x + (seq.int(0L, nx - 1L, by = stride) + 0.5) * spacing
y_values <- origin_y + (seq.int(0L, ny - 1L, by = stride) + 0.5) * spacing
grid <- CJ(x = x_values, y = y_values)
field[, r_active := r_active_small + r_active_large]
field <- field[, .(x, y, r_total, K_total, r_active, nutrient, time_hours)]
setkey(field, x, y)
setkey(grid, x, y)
plot_data <- field[grid]
for (column in c("r_total", "K_total", "r_active", "nutrient")) {
  set(plot_data, which(is.na(plot_data[[column]])), column, 0.0)
}
time_hours <- max(field$time_hours)
plot_data[, total := r_total + K_total]
plot_data[, r_share := fifelse(total > 1.0e-3, r_total / total, NA_real_)]

ink <- "#20252b"
quiet <- "#59616a"
blue <- "#2f6f9f"
orange <- "#d17a22"
panel_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(color = ink, face = "bold", size = 11),
    axis.title = element_text(color = quiet),
    axis.text = element_text(color = quiet),
    legend.title = element_text(color = quiet),
    legend.text = element_text(color = quiet)
  )
x_limits <- if (is.na(zoom_radius)) NULL else c(-zoom_radius, zoom_radius)
y_limits <- if (is.na(zoom_radius)) NULL else c(-zoom_radius, zoom_radius)

spatial_plot <- function(column, title, low, high, upper = 1.0) {
  ggplot(plot_data, aes(x = x, y = y, fill = .data[[column]])) +
    geom_raster() +
    scale_fill_gradient(
      low = low, high = high, limits = c(0.0, upper), oob = squish,
      name = "density"
    ) +
    coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    labs(title = title, x = "x (ABM voxel units)", y = "y (ABM voxel units)") +
    panel_theme
}

p_r <- spatial_plot("r_total", "r cell density", "#f7fbff", blue)
p_k <- spatial_plot("K_total", "K cell density", "#fff8f0", orange)
p_nutrient <- spatial_plot(
  "nutrient", "Effective nutrient supply", "#ffffff", "#30343b"
)
active_values <- plot_data[r_active > 0.0, r_active]
active_upper <- if (length(active_values)) {
  max(as.numeric(quantile(active_values, 0.995)), 1.0e-6)
} else {
  1.0
}
p_active <- spatial_plot(
  "r_active",
  sprintf("Activated r density (%gx migration)", activated_multiplier),
  "#f7fbff", blue,
  active_upper
)
p_share <- ggplot(plot_data, aes(x = x, y = y, fill = r_share)) +
  geom_raster() +
  scale_fill_gradient2(
    low = orange, mid = "#f4f1eb", high = blue, midpoint = 0.5,
    limits = c(0.0, 1.0), na.value = "white", name = "r share"
  ) +
  coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE) +
  labs(
    title = "Local r share of r + K",
    x = "x (ABM voxel units)", y = "y (ABM voxel units)"
  ) +
  panel_theme

plot_data[, radius_bin := floor(sqrt(x * x + y * y) / 8.0) * 8.0 + 4.0]
radial <- plot_data[, .(
  r = mean(r_total),
  K = mean(K_total)
), by = radius_bin]
radial <- radial[r + K > 1.0e-8]
radial_long <- melt(
  radial, id.vars = "radius_bin", variable.name = "population",
  value.name = "density"
)
p_radial <- ggplot(
  radial_long,
  aes(x = radius_bin, y = density, color = population, linetype = population)
) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c(r = blue, K = orange)) +
  scale_linetype_manual(values = c(r = "solid", K = "dashed")) +
  scale_x_continuous(expand = expansion(mult = c(0.0, 0.03))) +
  scale_y_continuous(limits = c(0.0, NA), expand = expansion(mult = c(0.0, 0.05))) +
  labs(
    title = "Azimuthal mean density", x = "Radius from lesion centre",
    y = "Mean density", color = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#d8dce1", linewidth = 0.3),
    plot.title = element_text(color = ink, face = "bold", size = 11),
    axis.title = element_text(color = quiet),
    axis.text = element_text(color = quiet),
    legend.position = "top"
  )

figure <- (p_r + p_k + p_nutrient) / (p_active + p_share + p_radial) +
  plot_annotation(
    title = sprintf(
      "ATCG3D structured PDE spatial distribution at %.0f hours", time_hours
    ),
    subtitle = paste(
      sprintf(
        "%d x %d simulation grid; displayed field sampled every %d site%s.",
        nx, ny, stride, if (stride == 1L) "" else "s"
      ),
      if (is.na(zoom_radius)) "" else sprintf(
        "Spatial panels show the central +/- %g voxel window.", zoom_radius
      ),
      "Density triggers activation; the stored clock controls return to baseline."
    ),
    theme = theme(
      plot.title = element_text(color = ink, face = "bold", size = 17),
      plot.subtitle = element_text(color = quiet, size = 10)
    )
  )

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
ggsave(output, figure, width = 16, height = 9, dpi = 180, bg = "white")

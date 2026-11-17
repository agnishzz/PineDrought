# =============================================================================
# Combine mortality–climate figures into one panel
#   Row 1: a = Pooled observed            b = Pooled simulated
#   Row 2: c = Spain observed (ET0)        d = Spain observed (CWB)   e = Spain simulated
#   Row 3: f = Switzerland observed       g = Switzerland simulated
#
# Spain observed is an additive model, so both terms are shown.
# fig_spain_obs.rds must contain BOTH panels: either the two-panel patchwork
# object, or list(p_eto, p_cwb). Element 1 = ET0, element 2 = CWB.
# =============================================================================

library(ggplot2)
library(patchwork)

# ---- 1. Panel definitions ---------------------------------------------------
# which: element to take from a multi-panel RDS | pos: "tr" or "tl" for stats text
base_dir <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/results_new1"

panels <- list(
  a = list(file = file.path(base_dir, "Pooled_ObsMor_F",    "fig_pooled_obs.rds"), type = "Observed",  var = "Prec", timing = "YL", window = "SS"),
  b = list(file = file.path(base_dir, "Pooled_SimMor",      "fig_pooled_sim.rds"), type = "Simulated", var = "Prec", timing = "YL", window = "SS"),
  c = list(file = file.path(base_dir, "Spain_ObsMor_F",     "fig_spain_obs.rds"),  type = "Observed",  var = "Eto",  timing = "AN", window = "SS", which = 1, pos = "tl"),
  d = list(file = file.path(base_dir, "Spain_ObsMor_F",     "fig_spain_obs.rds"),  type = "Observed",  var = "CWB",  timing = "AN", window = "SS", which = 2),
  e = list(file = file.path(base_dir, "Spain_SimMor",       "fig_spain_sim.rds"),  type = "Simulated", var = "Prec", timing = "AN", window = "SS"),
  f = list(file = file.path(base_dir, "Switzerland_ObsMor", "fig_swiss_obs.rds"),  type = "Observed",  var = "Prec", timing = "YL", window = "SS"),
  g = list(file = file.path(base_dir, "Switzerland_SimMor", "fig_swiss_sim.rds"),  type = "Simulated", var = "CWB",  timing = "YL", window = "SS")
)

# Stop early if any file is missing
missing <- vapply(panels, function(d) !file.exists(d$file), logical(1))
if (any(missing)) stop("Missing RDS file(s):\n", paste(unique(sapply(panels[missing], `[[`, "file")), collapse = "\n"))

# ---- 2. Full names for the x-axis (CHECK THESE) ----------------------------
var_names <- c(
  Prec = "Precipitation",
  CWB  = "Climatic water balance",
  Eto  = "Reference evapotranspiration"
)
timing_names <- c(
  YL = "previous year",   # <-- check
  AN = "current year"     # <-- check
)
window_names <- c(
  SS = "midsummer-early autumn"   # <-- check
)

make_xlab <- function(var, timing, window) {
  sprintf("%s anomaly\n(%s, %s)",
          var_names[[var]], timing_names[[timing]], window_names[[window]])
}

# ---- 3. Pull R2m and p from the existing subtitle ---------------------------
# Subtitle format: "b = -0.434, p = 0.014 | R2m = 0.10, R2c = 0.49"
get_stats <- function(p) {
  s <- p$labels$subtitle
  if (is.null(s) && exists("get_labs", asNamespace("ggplot2"))) s <- ggplot2::get_labs(p)$subtitle
  if (is.null(s)) stop("No subtitle with b/p/R2m found in one of the plots")
  p_str <- trimws(sub(".*p = ([^|]+?)\\s*\\|.*", "\\1", s))
  r2m   <- as.numeric(sub(".*R2m = ([0-9.]+).*", "\\1", s))
  list(p_str = p_str, r2m = r2m)
}

# p-value as plotmath: < 0.001, < 0.01, < 0.05, < 0.1, otherwise = value
p_plotmath <- function(p_str) {
  p <- as.numeric(sub("^<\\s*", "", p_str))
  if (p < 0.001)      "italic(p) < '0.001'"
  else if (p < 0.01)  "italic(p) < '0.01'"
  else if (p < 0.05)  "italic(p) < '0.05'"
  else if (p < 0.1)   "italic(p) < '0.1'"
  else                sprintf("italic(p) == '%.2f'", p)
}

# ---- 4. Load and restyle one panel ------------------------------------------
load_plot <- function(def) {
  obj <- readRDS(def$file)
  if (!is.null(def$which)) {
    obj <- obj[[def$which]]                 # pick one plot from a multi-plot file
  } else if (!inherits(obj, "ggplot")) {
    if (length(obj) == 1) {
      obj <- obj[[1]]                       # unwrap a one-element list
    } else {
      stop(basename(def$file), " contains ", length(obj), " plots: ",
           paste(names(obj), collapse = ", "),
           "\nAdd which = <number> to that panel's definition.")
    }
  }
  if (!inherits(obj, "ggplot")) stop(basename(def$file), " did not give a ggplot")
  obj
}

restyle <- function(def) {
  p  <- load_plot(def)
  st <- get_stats(p)
  
  r2_lab <- sprintf("italic(R)[m]^2 == '%.2f'", st$r2m)
  p_lab  <- p_plotmath(st$p_str)
  
  # Place the two lines in data units, as fractions of the panel's axis ranges
  pp <- ggplot_build(p)$layout$panel_params[[1]]
  xr <- pp$x.range; yr <- pp$y.range
  left <- identical(def$pos, "tl")
  x_at <- if (left) xr[1] + 0.03 * diff(xr) else xr[2] - 0.03 * diff(xr)
  hj   <- if (left) 0 else 1
  y1   <- yr[2] - (if (left) 0.14 else 0.02) * diff(yr)   # top-left: below the letter
  y2   <- y1 - 0.09 * diff(yr)
  
  p +
    labs(title = NULL, subtitle = NULL,
         x = make_xlab(def$var, def$timing, def$window),
         y = paste(def$type, "mortality (%)")) +
    annotate("text", x = x_at, y = y1, label = r2_lab, parse = TRUE,
             hjust = hj, vjust = 1, size = 3.8) +
    annotate("text", x = x_at, y = y2, label = p_lab, parse = TRUE,
             hjust = hj, vjust = 1, size = 3.8) +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      panel.grid       = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis lines so the panels still have a frame of reference without grid
      axis.line  = element_line(colour = "grey30", linewidth = 0.4),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.4)
    )
}

plots <- lapply(panels, restyle)

# ---- 5. Assemble: 2 panels in rows 1 and 3, 3 panels in row 2 --------------
design <- "
AAABBB
CCDDEE
FFFGGG
"

combined <- wrap_plots(plots, design = design) +
  # Panel letters, in the same order as `panels` (edit text freely)
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c-i)", "(c-ii)", "(d)", "(e)", "(f)"))) &
  theme(plot.tag          = element_text(face = "bold", size = 14),
        plot.tag.location = "panel",     # letter inside the plot area
        plot.tag.position = "topleft")

# ---- 6. Export --------------------------------------------------------------
ggsave(file.path(base_dir, "mortality_climate_panel3.png"), combined,
       width = 11, height = 11, dpi = 300, bg = "white")
# Vector version: points, lines and text stay editable in CorelDRAW
ggsave(file.path(base_dir, "mortality_climate_panel2.pdf"), combined,
       width = 12, height = 12, device = cairo_pdf, bg = "white")
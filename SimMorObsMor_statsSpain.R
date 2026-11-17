library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)

data <- readxl::read_xlsx("C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/NBR_Mor_Spain.xlsx")
head(data)

# =============================================================================
# SimMor (y) vs each remote-sensing-derived annual mortality series (x)
# - One linear regression per method: SimMor ~ ObsMor_method
# - R2, p-value, slope, intercept, RMSE, bias, Pearson r, n
# - Scatter + fitted lines coloured by method (+ faceted version)
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})

# ---- Input ------------------------------------------------------------------
# `data` = the combined table with Site, Year, A_*/B_* method columns and SimMor.
# If it is not already in your session, read it, e.g.:
# data <- read.csv("C:/.../spain_mortality_output/site_annual_mortpct_ALL_METHODS_long.csv")
out_dir <- "C:/Users/kumar_das/OneDrive - CZU v Praze/R_etal/PineDrought/simmor_vs_obsmor_output"
dir.create(out_dir, showWarnings = FALSE)

data <- as.data.frame(data)
y_col   <- "SimMor"
id_cols <- c("Site", "TotalTrees", "Dead", "PercentDead_census", "Year", "Plot")
x_cols  <- setdiff(names(data), c(id_cols, y_col))
x_cols  <- x_cols[sapply(data[x_cols], is.numeric)]
cat("Methods compared:", paste(x_cols, collapse = ", "), "\n\n")

# ---- Long format ------------------------------------------------------------
keep <- intersect(c("Site", "Plot", "Year"), names(data))
long <- do.call(rbind, lapply(x_cols, function(m)
  data.frame(data[keep], Method = m, ObsMor = data[[m]], SimMor = data[[y_col]])))
long <- long[complete.cases(long[, c("ObsMor", "SimMor")]), ]
long$Method <- factor(long$Method, levels = x_cols)
long$Scheme <- ifelse(startsWith(as.character(long$Method), "A_"), "A: plot signal",
                      ifelse(startsWith(as.character(long$Method), "B_"), "B: site-mean signal", "Other"))

# ---- Regression statistics per method ----------------------------------------
fit_stats <- function(d) {
  f <- lm(SimMor ~ ObsMor, data = d)
  s <- summary(f)
  cf <- coef(s)
  p  <- if (nrow(cf) > 1) cf[2, 4] else NA
  data.frame(
    n         = nrow(d),
    R2        = s$r.squared,
    adj_R2    = s$adj.r.squared,
    p_value   = p,
    slope     = unname(coef(f)[2]),
    intercept = unname(coef(f)[1]),
    pearson_r = cor(d$ObsMor, d$SimMor),
    spearman_rho = cor(d$ObsMor, d$SimMor, method = "spearman"),
    RMSE      = sqrt(mean((d$SimMor - d$ObsMor)^2)),   # vs 1:1 line
    bias      = mean(d$SimMor - d$ObsMor))             # SimMor - ObsMor
}

stats <- do.call(rbind, lapply(split(long, long$Method), function(d)
  cbind(Method = as.character(d$Method[1]), fit_stats(d))))
stats <- stats[order(-stats$R2), ]
stats$signif <- cut(stats$p_value, c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", ".", "ns"))
write.csv(stats, file.path(out_dir, "SimMor_vs_methods_stats.csv"), row.names = FALSE)

print_tbl <- stats
num <- sapply(print_tbl, is.numeric)
print_tbl[num] <- lapply(print_tbl[num], function(v) signif(v, 3))
print(print_tbl, row.names = FALSE)

# Same statistics per method within each site (optional, n = years per site)
if ("Site" %in% names(long)) {
  stats_site <- do.call(rbind, lapply(split(long, list(long$Method, long$Site), drop = TRUE), function(d)
    if (nrow(d) >= 3 && sd(d$ObsMor) > 0 && sd(d$SimMor) > 0)
      cbind(Method = as.character(d$Method[1]), Site = d$Site[1], fit_stats(d))))
  write.csv(stats_site, file.path(out_dir, "SimMor_vs_methods_stats_by_site.csv"), row.names = FALSE)
}

# ---- Colours: one per method; A = dark, B = light shade of the same hue -------
variants <- unique(sub("^[AB]_", "", x_cols))
hues <- setNames(seq(15, 375, length.out = length(variants) + 1)[seq_along(variants)], variants)
method_cols <- setNames(sapply(x_cols, function(m) {
  v <- sub("^[AB]_", "", m)
  if (startsWith(m, "B_")) hcl(hues[v], c = 55, l = 75) else hcl(hues[v], c = 85, l = 42)
}), x_cols)

# p-value text: "= 0.012" or "< 0.001"
p_txt <- function(p) ifelse(p < 0.001, "< 0.001", paste("=", formatC(p, format = "f", digits = 3)))

# Legend labels with R² and p (plotmath, so the superscript renders on any device)
lab <- lapply(x_cols, function(m) {
  i <- match(m, stats$Method)
  bquote(.(m) ~ "(" * R^2 == .(sprintf("%.2f", stats$R2[i])) * "," ~ p ~ .(p_txt(stats$p_value[i])) * ")")
})

lims <- range(c(long$ObsMor, long$SimMor), na.rm = TRUE)

# ---- Figure 1: all methods in one panel ---------------------------------------
p1 <- ggplot(long, aes(ObsMor, SimMor, colour = Method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, linewidth = 0.7) +
  scale_colour_manual(values = method_cols, breaks = x_cols, labels = lab) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = "Remote-sensing-derived annual mortality (%)",
       y = "Simulated annual mortality, SimMor (%)",
       colour = "Method", title = "SimMor vs remote-sensing mortality allocation methods",
       subtitle = "Dashed line = 1:1") +
  theme_bw() +
  theme(legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2.5)))
p1
ggsave(file.path(out_dir, "SimMor_vs_methods_allinone.png"), p1, width = 11, height = 7.5, dpi = 300)

# ---- Figure 2: one facet per method, R² and p annotated --------------------------
ann <- stats
ann$Method <- factor(ann$Method, levels = x_cols)
ann$label  <- sprintf("atop(R^2 == %.2f, p ~ '%s' ~~ n == %d)", ann$R2, p_txt(ann$p_value), ann$n)
shape_var <- if ("Site" %in% names(long)) "Site" else NULL

p2 <- ggplot(long, aes(ObsMor, SimMor, colour = Method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(shape = if (!is.null(shape_var)) .data[[shape_var]] else NULL), alpha = 0.75, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, linewidth = 0.7, fill = "grey80") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, parse = TRUE,
            hjust = -0.05, vjust = 1.1, size = 2.8) +
  facet_wrap(~ Method, ncol = 4) +
  scale_colour_manual(values = method_cols, guide = "none") +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = "Remote-sensing-derived annual mortality (%)",
       y = "Simulated annual mortality, SimMor (%)", shape = shape_var) +
  theme_bw() + theme(legend.position = "bottom", strip.text = element_text(size = 8))
p2
ggsave(file.path(out_dir, "SimMor_vs_methods_facets.png"), p2, width = 12, height = 10, dpi = 300)

# ---- Figure 3: R² bar chart per method -----------------------------------------
stats$Method <- factor(stats$Method, levels = rev(stats$Method))
p3 <- ggplot(stats, aes(R2, Method, fill = Method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(sprintf("%.2f", R2), " ", signif)), hjust = -0.1, size = 3) +
  scale_fill_manual(values = method_cols, guide = "none") +
  scale_x_continuous(limits = c(0, max(1, max(stats$R2) * 1.15)), expand = c(0, 0)) +
  labs(x = expression(R^2 ~ "(SimMor ~ method)"), y = NULL,
       caption = "*** p<0.001, ** p<0.01, * p<0.05, . p<0.1, ns p >= 0.1") +
  theme_bw()
ggsave(file.path(out_dir, "SimMor_vs_methods_R2.png"), p3, width = 7, height = 5, dpi = 300)

cat("\nOutputs written to", normalizePath(out_dir), "\n")

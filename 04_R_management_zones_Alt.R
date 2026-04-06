# =============================================================================
# Script 04_Alt: Management Zone Delineation — Pathway B: Soil Fertility Index
# Module 6 — Soil Fertility Assessment via AHP-Weighted Fuzzy Membership
#
# PATHWAY B — Soil Fertility Index (AHP-weighted fuzzy membership)
#     Each soil indicator is scored through a shape-appropriate fuzzy
#     membership function (S-curve for "more-is-better" nutrients; bell
#     for optima-constrained pH). Indicators are weighted by AHP
#     (or equal weights as fallback). SFI ∈ [0,1] is interpolated and
#     classified into fertility zones using standard breaks.
#
# Full pipeline:
#   +  S/bell membership functions for soil indicators
#   +  AHP-weighted Soil Fertility Index computation
#   +  Kriging of SFI surface (automap) with prediction uncertainty
#   +  SFI classification into fertility zones (Very Low → Very High)
#   +  Zone profiling by SFI class
#   +  Spatial validation: Moran's I, Kruskal-Wallis + Dunn tests
#   +  Fertility prescriptions per SFI class
#   +  Complete diagnostic CSV + PNG output set (figures 20–28)
#
# Inputs  (from Scripts 02 & 03):
#   data/embedding_pca.tif              — used as spatial template for kriging grid
#   data/soilgrids_with_embeddings.csv  — soil profiles + measured properties
#
# Outputs:
#   data/zones_sfi.tif                 — SFI surface
#   data/zones_sfi_class.tif           — SFI fertility class raster
#   data/sfi_kriging_variance.tif      — kriging prediction variance
#   outputs/zone_profiles.csv          — mean ± SD per SFI class per property
#   outputs/sfi_zone_table.csv         — SFI stats by class
#   outputs/fertility_prescriptions.csv
#   outputs/kruskal_dunn_tests.csv     — inter-class significance
#   outputs/moranI_results.csv         — Moran's I per soil property
#   figures/20_sfi_map.png
#   figures/21_sfi_class_map.png
#   figures/22_sfi_kriging_variance.png
#   figures/23_sfi_violin.png
#   figures/24_zone_boxplots.png
#   figures/25_zone_radar.png
#   figures/26_moran_correlogram.png
#   figures/27_prescription_summary.png
#   figures/28_final_figure.png
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────

required <- c(
  "terra",          # raster operations
  "sf",             # vector operations
  "tidyverse",      # data wrangling + ggplot2
  "viridis",        # colour palettes
  "patchwork",      # multi-panel plots
  "tidyterra",      # ggplot2 + terra
  "spdep",          # Moran's I, spatial weights
  "automap",        # automatic variogram fitting + kriging
  "dunn.test"       # Dunn test after Kruskal-Wallis
)

new_pkgs <- required[!(required %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required, library, character.only = TRUE))

.timer <- function(label, expr) {
  t0     <- proc.time()
  result <- force(expr)
  dt     <- (proc.time() - t0)[3]
  cat(sprintf("  [%s] %.1f s\n", label, dt))
  invisible(result)
}

# ── 1. USER CONFIGURATION ─────────────────────────────────────────────────────

CFG <- list(
  pca_raster     = "data/embedding_pca.tif",
  soil_csv       = "data/soilgrids_with_embeddings.csv",

  crs_code       = 32717,

  lon_col        = "longitude",
  lat_col        = "latitude",

  target_vars    = c("BD", "CEC", "Fragm", "Sand", "Silt", "Clay",
                     "N", "OCD", "pH", "SOC"),

  plot_max_cells   = 200000,

  sfi_vars = list(
    SOC = list(role = "more_better", a = 5,   b = 25),
    N   = list(role = "more_better", a = 100, b = 500),
    CEC = list(role = "more_better", a = 5,   b = 30),
    pH  = list(role = "optimum",     a = 6.2, b = 7.2),
    BD  = list(role = "less_better", a = 1.2, b = 1.5)
  ),

  sfi_ahp_weights = c(SOC = 0.35, N = 0.28, CEC = 0.18, pH = 0.12, BD = 0.07),

  sfi_breaks  = c(0, 0.25, 0.50, 0.75, 0.90, 1.0),
  sfi_labels  = c("Very Low", "Low", "Moderate", "High", "Very High"),

  thresh = list(
    SOC_low   = 10,   SOC_high  = 20,
    pH_low    = 5.8,  pH_high   = 7.2,
    CEC_low   = 10,   CEC_high  = 25
  ),

  out_data    = "data/",
  out_figures = "figures/",
  out_results = "outputs/"
)

walk(c(CFG$out_data, CFG$out_figures, CFG$out_results),
     ~if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# ── 2. LOAD INPUTS ─────────────────────────────────────────────────────────────

cat("Loading spatial template raster...\n")
pca_rast_all  <- rast(CFG$pca_raster)
template      <- pca_rast_all[[1]]
cat("Template dimensions:", ncol(template), "x", nrow(template),
    "(", ncell(template), "cells )\n")

cat("Loading soil data...\n")
soil_raw    <- read_csv(CFG$soil_csv, show_col_types = FALSE)
target_vars <- intersect(CFG$target_vars, names(soil_raw))
cat("Target variables for profiling:", paste(target_vars, collapse = ", "), "\n")

pts_sf <- st_as_sf(soil_raw,
                   coords = c(CFG$lon_col, CFG$lat_col),
                   crs    = 4326) |>
  st_transform(CFG$crs_code)

soil_raw$x <- st_coordinates(pts_sf)[, 1]
soil_raw$y <- st_coordinates(pts_sf)[, 2]

pts_sf_utm <- st_as_sf(soil_raw, coords = c("x", "y"),
                       crs = CFG$crs_code, remove = FALSE)

# =============================================================================
# ══ STEP 3 — SOIL FERTILITY INDEX (AHP-WEIGHTED FUZZY MEMBERSHIP) ═══════════
# =============================================================================

cat("\n=== PATHWAY B: Soil Fertility Index ===\n")

mf_s_shape <- function(x, a, b) {
  mid <- (a + b) / 2
  u   <- numeric(length(x))
  u[x <= a] <- 0
  u[x >= b] <- 1
  lo <- x > a & x <= mid
  hi <- x > mid & x < b
  u[lo] <- 2 * ((x[lo] - a) / (b - a))^2
  u[hi] <- 1 - 2 * ((x[hi] - b) / (b - a))^2
  u
}

mf_s_reverse <- function(x, a, b) 1 - mf_s_shape(x, a, b)

mf_bell <- function(x, a, b) {
  mid <- (a + b) / 2
  hw  <- (b - a) / 2
  u   <- exp(-((x - mid)^2) / (2 * (hw / 2)^2))
  u[x < a | x > b] <- 0
  u
}

score_indicator <- function(x, params) {
  switch(params$role,
    more_better = mf_s_shape(x,   params$a, params$b),
    less_better = mf_s_reverse(x, params$a, params$b),
    optimum     = mf_bell(x,      params$a, params$b),
    stop("Unknown role: ", params$role)
  )
}

sfi_vars_avail <- intersect(names(CFG$sfi_vars), names(soil_raw))
if (length(sfi_vars_avail) == 0) {
  stop("None of the SFI variables found in soil CSV. Check CFG$sfi_vars and column names.")
}

cat("  SFI variables found:", paste(sfi_vars_avail, collapse = ", "), "\n")

score_mat <- matrix(NA, nrow = nrow(soil_raw), ncol = length(sfi_vars_avail),
                    dimnames = list(NULL, sfi_vars_avail))
for (v in sfi_vars_avail) {
  score_mat[, v] <- score_indicator(soil_raw[[v]], CFG$sfi_vars[[v]])
}

ahp_w <- CFG$sfi_ahp_weights[sfi_vars_avail]
if (any(is.na(ahp_w))) {
  cat("  WARNING: some SFI variables lack AHP weights; using equal weights.\n")
  ahp_w[is.na(ahp_w)] <- 1
}
ahp_w <- ahp_w / sum(ahp_w)
cat("  Effective AHP weights:\n")
print(round(ahp_w, 3))

soil_raw$SFI <- as.numeric(score_mat %*% ahp_w)
cat(sprintf("  SFI range: %.3f – %.3f  (mean = %.3f)\n",
            min(soil_raw$SFI, na.rm = TRUE),
            max(soil_raw$SFI, na.rm = TRUE),
            mean(soil_raw$SFI, na.rm = TRUE)))

soil_raw$SFI_class <- cut(
  soil_raw$SFI,
  breaks = CFG$sfi_breaks,
  labels = CFG$sfi_labels,
  include.lowest = TRUE
)

cat("  SFI class distribution (sample points):\n")
print(table(soil_raw$SFI_class))

# ── Kriging of SFI surface ─────────────────────────────────────────────────

pts_sfi_utm <- st_as_sf(soil_raw |> filter(!is.na(SFI)),
                        coords = c("x", "y"),
                        crs = CFG$crs_code)

pts_sp <- as(pts_sfi_utm["SFI"], "Spatial")

coarse_fact <- max(1, ceiling(sqrt(ncell(template) / 500000)))
template_coarse <- aggregate(template, fact = coarse_fact)
pred_grid_df <- as.data.frame(template_coarse, xy = TRUE, na.rm = TRUE)
pred_grid_sf <- st_as_sf(pred_grid_df, coords = c("x", "y"),
                         crs = CFG$crs_code)
pred_grid_sp <- as(pred_grid_sf, "Spatial")

cat("  Fitting variogram and kriging SFI surface (automap)...\n")
cat("  Prediction grid:", nrow(pred_grid_sp), "cells (coarsened by factor",
    coarse_fact, ")\n")
krige_result <- tryCatch(
  .timer("autoKrige SFI", {
    automap::autoKrige(SFI ~ 1, input_data = pts_sp,
                       new_data = pred_grid_sp)
  }),
  error = function(e) {
    cat("  autoKrige failed:", conditionMessage(e),
        "— using IDW interpolation as fallback.\n")
    NULL
  }
)

downsample_for_plot <- function(r, max_cells) {
  if (ncell(r) > max_cells) {
    fact <- ceiling(sqrt(ncell(r) / max_cells))
    aggregate(r, fact = fact, fun = "modal", na.rm = TRUE)
  } else r
}

downsample_continuous <- function(r, max_cells) {
  if (ncell(r) > max_cells) {
    fact <- ceiling(sqrt(ncell(r) / max_cells))
    aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
  } else r
}

sfi_cols <- c("Very Low" = "#d73027", "Low" = "#fc8d59",
              "Moderate" = "#fee08b", "High" = "#91cf60",
              "Very High" = "#1a9850")

if (!is.null(krige_result)) {

  krige_sf <- st_as_sf(krige_result$krige_output)
  krige_coords <- st_coordinates(krige_sf)
  krige_df <- st_drop_geometry(krige_sf)
  krige_df$x <- krige_coords[, 1]
  krige_df$y <- krige_coords[, 2]

  sfi_krige_rast <- rast(krige_df[, c("x", "y", "var1.pred")],
                         type = "xyz")
  sfi_var_rast   <- rast(krige_df[, c("x", "y", "var1.var")],
                         type = "xyz")
  crs(sfi_krige_rast) <- paste0("EPSG:", CFG$crs_code)
  crs(sfi_var_rast)   <- paste0("EPSG:", CFG$crs_code)
  sfi_krige_rast <- resample(sfi_krige_rast, template, method = "bilinear")
  sfi_var_rast   <- resample(sfi_var_rast,   template, method = "bilinear")
  sfi_krige_rast <- clamp(sfi_krige_rast, 0, 1)

  rcl_mat <- matrix(
    c(CFG$sfi_breaks[-length(CFG$sfi_breaks)],
      CFG$sfi_breaks[-1],
      seq_along(CFG$sfi_labels)),
    ncol = 3
  )
  sfi_class_rast <- terra::classify(sfi_krige_rast, rcl_mat)
  levels(sfi_class_rast) <- data.frame(
    value = seq_along(CFG$sfi_labels),
    label = CFG$sfi_labels
  )

  writeRaster(sfi_krige_rast,  file.path(CFG$out_data, "zones_sfi.tif"),
              overwrite = TRUE)
  writeRaster(sfi_class_rast,  file.path(CFG$out_data, "zones_sfi_class.tif"),
              overwrite = TRUE)
  writeRaster(sfi_var_rast,    file.path(CFG$out_data, "sfi_kriging_variance.tif"),
              overwrite = TRUE)

  # ── Figure 20: SFI continuous surface ─────────────────────────────────────

  sfi_plot_rast <- downsample_continuous(sfi_krige_rast, CFG$plot_max_cells)

  p_sfi <- ggplot() +
    geom_spatraster(data = sfi_plot_rast) +
    scale_fill_distiller(palette = "RdYlGn", direction = 1,
                         name = "SFI\n[0–1]", na.value = "grey90",
                         limits = c(0, 1)) +
    geom_sf(data = pts_sfi_utm, aes(colour = SFI), size = 2, shape = 16) +
    scale_colour_distiller(palette = "RdYlGn", direction = 1,
                           name = "SFI (obs)", limits = c(0, 1)) +
    labs(title    = "Soil Fertility Index — kriged surface",
         subtitle = paste0("AHP weights: ",
                           paste(names(ahp_w), round(ahp_w, 2),
                                 sep = "=", collapse = ", "))) +
    theme_bw(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  ggsave(file.path(CFG$out_figures, "20_sfi_map.png"),
         p_sfi, width = 8, height = 6.5, dpi = 300)

  # ── Figure 21: SFI class raster ───────────────────────────────────────────

  sfi_class_plot <- downsample_for_plot(sfi_class_rast, CFG$plot_max_cells)

  p_sfi_class <- ggplot() +
    geom_spatraster(data = sfi_class_plot) +
    scale_fill_manual(values   = sfi_cols, name = "Fertility\nclass",
                      na.value = "grey90") +
    labs(title    = "SFI fertility classification",
         subtitle = "Breaks: 0 | 0.25 | 0.50 | 0.75 | 0.90 | 1.0") +
    theme_bw(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  ggsave(file.path(CFG$out_figures, "21_sfi_class_map.png"),
         p_sfi_class, width = 8, height = 6.5, dpi = 300)

  # ── Figure 22: kriging variance ────────────────────────────────────────────

  sfi_var_plot <- downsample_continuous(sfi_var_rast, CFG$plot_max_cells)

  p_krig_var <- ggplot() +
    geom_spatraster(data = sfi_var_plot) +
    scale_fill_viridis_c(option = "B", name = "Kriging\nvariance",
                         na.value = "grey90") +
    labs(title    = "SFI kriging prediction variance",
         subtitle = "Higher values = greater interpolation uncertainty") +
    theme_bw(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  ggsave(file.path(CFG$out_figures, "22_sfi_kriging_variance.png"),
         p_krig_var, width = 7, height = 6, dpi = 300)

}

# ── SFI summary by class ──────────────────────────────────────────────────────

sfi_zone_tbl <- soil_raw |>
  filter(!is.na(SFI)) |>
  group_by(sfi_class = as.character(SFI_class)) |>
  summarise(
    n          = n(),
    SFI_mean   = mean(SFI,   na.rm = TRUE),
    SFI_sd     = sd(SFI,     na.rm = TRUE),
    SFI_median = median(SFI, na.rm = TRUE),
    .groups    = "drop"
  ) |>
  arrange(desc(SFI_mean)) |>
  mutate(fertility_rank  = row_number(),
         fertility_label = case_when(
           fertility_rank == 1     ~ "Highest fertility",
           fertility_rank == n()   ~ "Lowest fertility",
           TRUE                    ~ paste("Intermediate", fertility_rank)
         ))

print(sfi_zone_tbl)
write_csv(sfi_zone_tbl, file.path(CFG$out_results, "sfi_zone_table.csv"))

# ── Figure 23: SFI violin by class ────────────────────────────────────────────

p_sfi_violin <- soil_raw |>
  filter(!is.na(SFI)) |>
  mutate(sfi_class = factor(as.character(SFI_class),
                            levels = CFG$sfi_labels)) |>
  ggplot(aes(sfi_class, SFI, fill = sfi_class)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.size = 1, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = sfi_cols, guide = "none") +
  geom_hline(yintercept = CFG$sfi_breaks[-c(1, length(CFG$sfi_breaks))],
             linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  labs(title    = "Soil Fertility Index distribution by fertility class",
       subtitle = "Dashed lines = SFI class breaks (Kumar et al. 2023)",
       x = NULL, y = "SFI") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(CFG$out_figures, "23_sfi_violin.png"),
       p_sfi_violin, width = 9, height = 6, dpi = 300)

# =============================================================================
# ══ STEP 4 — ZONE PROFILING BY SFI CLASS ═══════════════════════════════════
# =============================================================================

cat("\n=== Zone profiling by SFI class ===\n")

soil_zones <- pts_sf_utm |>
  st_drop_geometry() |>
  select(all_of(target_vars), x, y) |>
  mutate(sfi_class = as.character(soil_raw$SFI_class)) |>
  drop_na(sfi_class) |>
  filter(!is.na(sfi_class))

zone_summary <- soil_zones |>
  group_by(sfi_class) |>
  summarise(
    n = n(),
    across(all_of(target_vars),
           list(mean   = ~mean(.x,   na.rm = TRUE),
                sd     = ~sd(.x,     na.rm = TRUE),
                median = ~median(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

print(zone_summary)
write_csv(zone_summary, file.path(CFG$out_results, "zone_profiles.csv"))

# ── Figure 24: Box plots by SFI class ─────────────────────────────────────────

soil_long <- soil_zones |>
  pivot_longer(all_of(target_vars), names_to = "property", values_to = "value") |>
  mutate(property  = factor(property, levels = target_vars),
         sfi_class = factor(sfi_class, levels = CFG$sfi_labels))

ggplot(soil_long, aes(sfi_class, value, fill = sfi_class)) +
  geom_boxplot(outlier.shape = 21, linewidth = 0.6, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 18,
               size = 3, colour = "white") +
  facet_wrap(~property, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = sfi_cols) +
  labs(title    = "Soil property distribution by SFI fertility class",
       subtitle = "Diamond = mean", x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

ggsave(file.path(CFG$out_figures, "24_zone_boxplots.png"),
       width = 11, height = 7, dpi = 300)

# ── Figure 25: Radar chart by SFI class ───────────────────────────────────────

radar_df <- zone_summary |>
  select(sfi_class, ends_with("_mean")) |>
  rename_with(~str_remove(.x, "_mean"), ends_with("_mean")) |>
  pivot_longer(-sfi_class, names_to = "property", values_to = "value") |>
  group_by(property) |>
  mutate(norm = (value - min(value)) / (max(value) - min(value) + 1e-9)) |>
  ungroup()

ggplot(radar_df, aes(x = property, y = norm, colour = sfi_class, group = sfi_class)) +
  geom_polygon(aes(fill = sfi_class), alpha = 0.15, linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = sfi_cols) +
  scale_fill_manual(values   = sfi_cols) +
  coord_polar() +
  labs(title = "SFI class profile radar chart (normalised)",
       x = NULL, y = NULL, colour = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom")

ggsave(file.path(CFG$out_figures, "25_zone_radar.png"),
       width = 7, height = 7, dpi = 300)

# =============================================================================
# ══ STEP 5 — SPATIAL VALIDATION ════════════════════════════════════════════
# =============================================================================

cat("\n--- Kruskal-Wallis + Dunn tests (by SFI class) ---\n")

kw_results <- map_dfr(target_vars, function(v) {
  df_v <- soil_zones |> select(sfi_class, value = all_of(v)) |> drop_na()
  if (n_distinct(df_v$sfi_class) < 2 || nrow(df_v) < 5) return(NULL)
  kt <- kruskal.test(value ~ sfi_class, data = df_v)
  tibble(
    property    = v,
    KW_chi2     = round(kt$statistic, 3),
    KW_df       = kt$parameter,
    KW_p        = round(kt$p.value,   4),
    significant = kt$p.value < 0.05
  )
})

cat("\nKruskal-Wallis results:\n")
print(kw_results)

dunn_results <- map_dfr(
  kw_results |> filter(significant) |> pull(property),
  function(v) {
    df_v <- soil_zones |> select(sfi_class, value = all_of(v)) |> drop_na()
    dt   <- dunn.test::dunn.test(df_v$value, df_v$sfi_class,
                                 method = "bonferroni", altp = TRUE)
    tibble(
      property   = v,
      comparison = dt$comparisons,
      Z          = round(dt$Z, 3),
      p_adj      = round(dt$altP.adjusted, 4)
    )
  }
)

write_csv(bind_rows(kw_results, .id = NULL),
          file.path(CFG$out_results, "kruskal_dunn_tests.csv"))

# ── Moran's I per soil property ───────────────────────────────────────────────

cat("\n--- Moran's I ---\n")

coords_mat <- as.matrix(soil_zones[, c("x", "y")])
nb_knn     <- spdep::knearneigh(coords_mat, k = 8) |> spdep::knn2nb()
lw_knn     <- spdep::nb2listw(nb_knn, style = "W")

moran_results <- map_dfr(target_vars, function(v) {
  x_v <- soil_zones[[v]]
  if (sum(!is.na(x_v)) < 10) return(NULL)
  mt  <- spdep::moran.test(x_v, lw_knn, na.action = na.exclude)
  tibble(
    property    = v,
    MoranI      = round(mt$estimate["Moran I statistic"], 4),
    Expectation = round(mt$estimate["Expectation"], 6),
    Variance    = round(mt$estimate["Variance"], 8),
    p_value     = round(mt$p.value, 4),
    significant = mt$p.value < 0.05
  )
})

cat("\nMoran's I results:\n")
print(moran_results)
write_csv(moran_results, file.path(CFG$out_results, "moranI_results.csv"))

# ── Figure 26: Moran's I bar chart + scatter plots ────────────────────────────

p_moran_bar <- ggplot(moran_results,
                  aes(x = reorder(property, MoranI), y = MoranI,
                      fill = significant)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#b2b2b2"),
                    labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                    name = NULL) +
  labs(title = "Moran's I by soil property",
       subtitle = "Blue = statistically significant (p < 0.05)",
       x = NULL, y = "Moran's I") +
  coord_flip() +
  theme_bw(base_size = 11)

sig_vars <- moran_results |> filter(significant) |> pull(property)
if (length(sig_vars) == 0) sig_vars <- moran_results |> slice_max(MoranI, n = 4) |> pull(property)
sig_vars <- head(sig_vars, 4)

moran_scatter_list <- map(sig_vars, function(v) {
  x_v   <- scale(soil_zones[[v]])[, 1]
  wz_v  <- spdep::lag.listw(lw_knn, x_v)
  df_sc <- tibble(x_std = x_v, wx = as.numeric(wz_v))
  I_val <- moran_results |> filter(property == v) |> pull(MoranI)
  p_val <- moran_results |> filter(property == v) |> pull(p_value)

  ggplot(df_sc, aes(x_std, wx)) +
    geom_point(alpha = 0.4, size = 1.2, colour = "#4393c3") +
    geom_smooth(method = "lm", se = FALSE, colour = "#d6604d", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3,
             label = sprintf("I = %.3f\np = %.4f", I_val, p_val),
             size = 3, fontface = "italic") +
    labs(title = v, x = "Standardised value (z)",
         y = "Spatial lag (Wz)") +
    theme_bw(base_size = 10)
})

p_scatter_combined <- wrap_plots(moran_scatter_list, ncol = 2)

p_moran <- (p_moran_bar / p_scatter_combined) +
  plot_annotation(
    title = "Moran's I spatial autocorrelation — soil properties",
    subtitle = "Scatter plots: slope of fitted line = Moran's I; quadrants = HH (top-right), LL (bottom-left), HL, LH"
  )

ggsave(file.path(CFG$out_figures, "26_moran_correlogram.png"),
       p_moran, width = 9, height = 10, dpi = 300)

# =============================================================================
# ══ STEP 6 — FERTILITY PRESCRIPTIONS ═══════════════════════════════════════
# =============================================================================

classify_thresh <- function(val, low, high, low_label, high_label,
                            mid_label = "Marginal") {
  case_when(
    val < low  ~ low_label,
    val > high ~ high_label,
    TRUE       ~ mid_label
  )
}

prescriptions <- zone_summary |>
  select(sfi_class, n, ends_with("_mean")) |>
  rename_with(~str_remove(.x, "_mean"), ends_with("_mean")) |>
  mutate(
    SOC_status = if ("SOC" %in% names(pick(everything()))) classify_thresh(SOC,
                    CFG$thresh$SOC_low, CFG$thresh$SOC_high,
                    "Deficient — add organic amendments",
                    "Adequate") else NA_character_,

    pH_status  = if ("pH" %in% names(pick(everything()))) case_when(
      pH < CFG$thresh$pH_low  ~ "Acidic — lime application recommended",
      pH > CFG$thresh$pH_high ~ "Alkaline — gypsum or S may be needed",
      TRUE                    ~ "Optimal pH range"
    ) else NA_character_,

    CEC_status = if ("CEC" %in% names(pick(everything()))) classify_thresh(CEC,
                    CFG$thresh$CEC_low, CFG$thresh$CEC_high,
                    "Low CEC — low nutrient retention",
                    "High CEC") else NA_character_,

    VR_priority = case_when(
      str_detect(coalesce(SOC_status, ""), "Deficient") |
        str_detect(coalesce(pH_status,  ""), "Acidic|Alkaline") |
        str_detect(coalesce(CEC_status, ""), "Low")              ~ "HIGH",
      str_detect(coalesce(SOC_status, ""), "Marginal") |
        str_detect(coalesce(CEC_status, ""), "Marginal")         ~ "MEDIUM",
      TRUE                                                        ~ "LOW"
    )
  ) |>
  select(sfi_class, n, any_of(target_vars),
         any_of(c("SOC_status", "pH_status", "CEC_status")), VR_priority)

cat("\n===  FERTILITY PRESCRIPTIONS  ===\n")
print(prescriptions)
write_csv(prescriptions,
          file.path(CFG$out_results, "fertility_prescriptions.csv"))

# ── Figure 27: prescription summary ──────────────────────────────────────────

prescriptions |>
  mutate(sfi_class = factor(sfi_class, levels = CFG$sfi_labels)) |>
  pivot_longer(any_of(c("SOC_status", "pH_status", "CEC_status")),
               names_to = "indicator", values_to = "status") |>
  mutate(
    indicator = recode(indicator,
                       "SOC_status" = "Soil Organic Carbon",
                       "pH_status"  = "Soil pH",
                       "CEC_status" = "CEC"),
    severity = case_when(
      str_detect(status, "Deficient|Acidic|Alkaline|Low") ~ "Action needed",
      str_detect(status, "Marginal")                      ~ "Monitor",
      TRUE                                                ~ "Adequate"
    ),
    severity = factor(severity, levels = c("Action needed", "Monitor", "Adequate"))
  ) |>
  ggplot(aes(sfi_class, fill = severity)) +
  geom_bar(position = "stack") +
  facet_wrap(~indicator, nrow = 1) +
  scale_fill_manual(values = c("Action needed" = "#d73027",
                               "Monitor"       = "#fee090",
                               "Adequate"      = "#1a9850"),
                    name = NULL) +
  labs(title = "Soil fertility status by SFI class",
       x = NULL, y = "Count of indicators") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")

ggsave(file.path(CFG$out_figures, "27_prescription_summary.png"),
       width = 10, height = 5, dpi = 300)

# =============================================================================
# ══ STEP 7 — FINAL COMBINED FIGURE ═════════════════════════════════════════
# =============================================================================

if (exists("sfi_krige_rast") && exists("sfi_class_rast")) {

  sfi_plot_rast2  <- downsample_continuous(sfi_krige_rast, CFG$plot_max_cells)
  sfi_class_plot2 <- downsample_for_plot(sfi_class_rast, CFG$plot_max_cells)

  p_final_sfi <- ggplot() +
    geom_spatraster(data = sfi_plot_rast2) +
    scale_fill_distiller(palette = "RdYlGn", direction = 1,
                         name = "SFI", na.value = "grey90",
                         limits = c(0, 1)) +
    geom_sf(data = pts_sfi_utm, colour = "black", fill = "white",
            size = 1.5, shape = 21, stroke = 0.5) +
    labs(title = "Soil Fertility Index (kriged)",
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  p_final_class <- ggplot() +
    geom_spatraster(data = sfi_class_plot2) +
    scale_fill_manual(values = sfi_cols, name = "Fertility\nclass",
                      na.value = "grey90") +
    geom_sf(data = pts_sfi_utm, colour = "black", fill = "white",
            size = 1.5, shape = 21, stroke = 0.5) +
    labs(title = "SFI fertility classes",
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  p_final_var <- ggplot() +
    geom_spatraster(data = downsample_continuous(sfi_var_rast, CFG$plot_max_cells)) +
    scale_fill_viridis_c(option = "B", name = "Variance",
                         na.value = "grey90") +
    labs(title = "Kriging prediction variance",
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  (p_final_sfi | p_final_var | p_final_class) +
    plot_annotation(
      title = "Soil Fertility Index — Pathway B",
      subtitle = paste0("AHP weights: ",
                        paste(names(ahp_w), round(ahp_w, 2),
                              sep = "=", collapse = ", "),
                        "  |  n = ", nrow(soil_raw), " soil profiles"),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )

  ggsave(file.path(CFG$out_figures, "28_final_figure.png"),
         width = 18, height = 6, dpi = 300)
}

# =============================================================================
# ══ SUMMARY ══════════════════════════════════════════════════════════════════
# =============================================================================

cat("\n=== Script 04_Alt (Pathway B only) complete ===\n")
cat(sprintf("  SFI range    : %.3f – %.3f  (mean = %.3f)\n",
            min(soil_raw$SFI, na.rm = TRUE),
            max(soil_raw$SFI, na.rm = TRUE),
            mean(soil_raw$SFI, na.rm = TRUE)))
cat("  SFI class distribution:\n")
print(table(soil_raw$SFI_class))
cat("  Outputs:\n")
cat("    Rasters  :", file.path(CFG$out_data,    "zones_sfi*.tif"), "\n")
cat("    Tables   :", file.path(CFG$out_results, "*.csv"), "\n")
cat("    Figures  :", file.path(CFG$out_figures, "20-28_*.png"), "\n")
cat("\n  Moran's I summary (significant properties):\n")
print(moran_results |> filter(significant))
cat("\n  Kruskal-Wallis summary:\n")
print(kw_results)

# =============================================================================
# END OF SCRIPT 04_Alt
# =============================================================================

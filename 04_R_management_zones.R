# =============================================================================
# Script 04: Management Zone Delineation from Satellite Embedding PCs
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────

required <- c(
  "terra",          # raster operations
  "sf",             # vector operations
  "tidyverse",      # data wrangling + ggplot2
  "ppclust",        # fuzzy c-means (fcm)
  "cluster",        # silhouette index
  "ClustGeo",       # spatially constrained hierarchical clustering
  "viridis",        # colour palettes
  "patchwork",      # multi-panel plots
  "tidyterra",      # ggplot2 + terra
  "factoextra",     # cluster diagnostics
  # --- NEW optimisation packages ---
  "ClusterR",       # mini-batch k-means (C++ backend)
  "matrixStats",    # fast row-wise operations (rowMaxs, etc.)
  "parallelDist"    # multi-threaded distance matrices (C++)
)

new_pkgs <- required[!(required %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required, library, character.only = TRUE))

# Helper: wall-clock timer
.timer <- function(label, expr) {
  t0 <- proc.time()
  result <- force(expr)
  dt <- (proc.time() - t0)[3]
  cat(sprintf("  [%s] %.1f s\n", label, dt))
  invisible(result)
}

# ── 1. USER CONFIGURATION ─────────────────────────────────────────────────────

CFG <- list(
  # Input files
  pca_raster     = "data/embedding_pca.tif",
  soil_csv       = "data/soilgrids_with_embeddings.csv",
  pred_raster    = "data/soil_predictions.tif",

  # CRS (EPSG:32717 — UTM Zone 17S, WGS84)
  crs_code       = 32717,

  # Coordinate columns (WGS84 lon/lat — will be reprojected to UTM)
  lon_col        = "longitude",
  lat_col        = "latitude",

  # Soil properties available in the CSV
  target_vars    = c("BD", "CEC", "Fragm", "Sand", "Silt", "Clay",
                     "N", "OCD", "pH", "SOC"),

  # Clustering settings
  n_pcs_cluster  = 10,
  k_min          = 2,
  k_max          = 10,
  k_final        = NULL,     # NULL = auto-select by silhouette peak

  # Fuzzy c-means parameters
  fuzziness      = 2.0,
  fcm_nstart     = 20,

  # ClustGeo parameters
  alpha_grid     = seq(0, 1, 0.1),
  geo_unit       = 1000,

  # Subsampling sizes  (OPTIMISED: added separate FCM subsample)
  sample_size_ksel  = 5000,   # k-selection elbow/silhouette
  sample_size_fcm   = 30000,  # FCM fitting (new — avoids 4M-row FCM)
  sample_size_sil   = 2000,   # silhouette evaluation (separate, smaller)
  sample_size_geo   = 3000,   # ClustGeo

  # Plot downsampling — max pixels rendered in ggplot maps
  plot_max_cells    = 200000,

  # Fertility thresholds
  thresh = list(
    SOC_low   = 10,
    SOC_high  = 20,
    pH_low    = 5.8,
    pH_high   = 7.2,
    CEC_low   = 10,
    CEC_high  = 25
  ),

  # Output directories
  out_data    = "data/",
  out_figures = "figures/",
  out_results = "outputs/"
)

walk(c(CFG$out_data, CFG$out_figures, CFG$out_results),
     ~if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# ── 2. LOAD INPUTS ─────────────────────────────────────────────────────────────

cat("Loading PC raster...\n")
pca_rast_all   <- rast(CFG$pca_raster)
n_pcs_cluster  <- min(CFG$n_pcs_cluster, nlyr(pca_rast_all))
pca_rast       <- pca_rast_all[[1:n_pcs_cluster]]
names(pca_rast) <- paste0("PC", 1:n_pcs_cluster)
cat("PCs used for clustering:", n_pcs_cluster,
    "(requested:", CFG$n_pcs_cluster, ")\n")

# Extract all non-NA pixel values as a matrix
# OPT: terra::as.data.frame with xy=TRUE is already efficient; keep as-is
df_all  <- .timer("Extract raster to data.frame", {
  as.data.frame(pca_rast, xy = TRUE, na.rm = TRUE)
})
xy_all  <- df_all[, c("x", "y")]
mat_all <- as.matrix(df_all[, paste0("PC", 1:n_pcs_cluster)])

# OPT: Pre-compute scale parameters once, apply to full matrix
center_vals <- colMeans(mat_all)
scale_vals  <- apply(mat_all, 2, sd)
mat_scl     <- scale(mat_all, center = center_vals, scale = scale_vals)

n_pixels <- nrow(mat_all)
cat("Pixels available for clustering:", n_pixels, "\n")

# Soil profile data
soil_raw <- read_csv(CFG$soil_csv, show_col_types = FALSE)
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

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — SELECT OPTIMAL NUMBER OF ZONES
# ─────────────────────────────────────────────────────────────────────────────

set.seed(99)
idx_sel <- sample(n_pixels, min(CFG$sample_size_ksel, n_pixels))
sub_scl <- mat_scl[idx_sel, ]

# Smaller subsample for silhouette (dist() on 5000 points = 25M entries = slow)
idx_sil <- sample(length(idx_sel), min(CFG$sample_size_sil, length(idx_sel)))
sub_sil <- sub_scl[idx_sil, ]

cat("Computing elbow and silhouette for k =",
    CFG$k_min, "to", CFG$k_max, "...\n")

k_range <- CFG$k_min:CFG$k_max

# OPT: Pre-compute distance matrix for silhouette using parallelDist (C++)
cat("  Computing distance matrix for silhouette (n =",
    nrow(sub_sil), ")...\n")
dist_sil <- .timer("parallelDist for silhouette", {
  parallelDist(sub_sil, method = "euclidean")
})

wss_vec <- sil_vec <- numeric(length(k_range))

.timer("k-selection loop", {
  for (i in seq_along(k_range)) {
    k <- k_range[i]

    # OPT: MiniBatchKmeans on the larger subsample for WSS
    mb <- MiniBatchKmeans(sub_scl, clusters = k, batch_size = 500,
                          num_init = 10, max_iters = 200,
                          init_fraction = 0.3, early_stop_iter = 10,
                          verbose = FALSE)
    # WSS from the mini-batch centroids — recompute on full subsample
    cl_full <- predict_MBatchKMeans(sub_scl, mb$centroids)
    wss_vec[i] <- sum(vapply(seq_len(k), function(j) {
      idx_j <- which(cl_full == j)
      if (length(idx_j) < 2) return(0)
      sum(sweep(sub_scl[idx_j, , drop = FALSE],
                2, mb$centroids[j, ])^2)
    }, numeric(1)))

    # Silhouette on the smaller subsample (fast)
    cl_sil <- predict_MBatchKMeans(sub_sil, mb$centroids)
    ss <- silhouette(cl_sil, dist_sil)
    sil_vec[i] <- mean(ss[, 3])
  }
})

sel_df <- tibble(k = k_range, WSS = wss_vec, Silhouette = sil_vec)

p_wss <- ggplot(sel_df, aes(k, WSS)) +
  geom_line(colour = "#1a7abf", linewidth = 1) +
  geom_point(size = 3, colour = "#1a7abf") +
  scale_x_continuous(breaks = k_range) +
  labs(title = "Elbow criterion",
       x = "Number of zones (k)", y = "Total within-cluster SS") +
  theme_bw(base_size = 11)

p_sil <- ggplot(sel_df, aes(k, Silhouette)) +
  geom_line(colour = "#e36b2d", linewidth = 1) +
  geom_point(size = 3, colour = "#e36b2d") +
  scale_x_continuous(breaks = k_range) +
  labs(title = "Average silhouette width",
       x = "Number of zones (k)", y = "Silhouette") +
  theme_bw(base_size = 11)

p_zone_sel <- (p_wss | p_sil) +
  plot_annotation(title = "Zone number selection — embedding PCs")
ggsave(file.path(CFG$out_figures, "11_zone_selection.png"),
       p_zone_sel, width = 11, height = 4.5, dpi = 300)

k_opt <- if (!is.null(CFG$k_final)) {
  CFG$k_final
} else {
  k_range[which.max(sil_vec)]
}
cat("Selected k =", k_opt, "(silhouette peak)\n")
cat(">>> Inspect figures/11_zone_selection.png and override CFG$k_final if needed.\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — FUZZY C-MEANS CLUSTERING
# ─────────────────────────────────────────────────────────────────────────────

cat("\nRunning fuzzy c-means (k =", k_opt, ", m =", CFG$fuzziness, ")...\n")
set.seed(42)

n_fcm <- min(CFG$sample_size_fcm, n_pixels)
if (n_fcm < n_pixels) {
  idx_fcm <- sample(n_pixels, n_fcm)
  mat_fcm <- mat_scl[idx_fcm, ]
  cat("  FCM fitted on subsample:", n_fcm, "of", n_pixels, "pixels\n")
} else {
  mat_fcm <- mat_scl
  cat("  FCM fitted on all", n_pixels, "pixels\n")
}

fcm_res <- .timer("FCM fitting", {
  fcm(
    x        = mat_fcm,
    centers  = k_opt,
    m        = CFG$fuzziness,
    nstart   = CFG$fcm_nstart,
    iter.max = 200
  )
})

cat("  Assigning full grid (", n_pixels, "pixels) to FCM centroids...\n")
centroids <- fcm_res$v   # k_opt x n_pcs matrix

.timer("Full-grid FCM assignment", {
  # Squared Euclidean distance from each pixel to each centroid
  # Result: n_pixels x k_opt matrix
  dist2_mat <- matrix(0, nrow = n_pixels, ncol = k_opt)
  for (j in seq_len(k_opt)) {
    diff_mat <- sweep(mat_scl, 2, centroids[j, ])
    dist2_mat[, j] <- rowSums(diff_mat^2)
  }

  # Fuzzy membership: u_ij = 1 / sum_l( (d2_ij / d2_il)^(1/(m-1)) )
  exp_val <- 1 / (CFG$fuzziness - 1)
  membership_mat <- matrix(0, nrow = n_pixels, ncol = k_opt)
  for (j in seq_len(k_opt)) {
    # Ratio d2_ij / d2_il for each l, summed
    ratio_sum <- rowSums(
      (dist2_mat[, j] / pmax(dist2_mat, .Machine$double.eps))^exp_val
    )
    membership_mat[, j] <- 1 / ratio_sum
  }
})

# Use matrixStats for row-wise max (C-level, no R loop)
max_member  <- rowMaxs(membership_mat)
hard_zones  <- max.col(membership_mat, ties.method = "first")
uncertainty <- 1 - max_member

df_all$zone_fcm    <- hard_zones
df_all$max_member  <- max_member
df_all$uncertainty <- uncertainty

member_col_names <- paste0("member_zone", 1:k_opt)
for (j in seq_len(k_opt)) {
  df_all[[member_col_names[j]]] <- membership_mat[, j]
}

# ── Convert to rasters (single template copy, direct cell assignment) ───

cat("  Building output rasters...\n")
template <- pca_rast[[1]]

# OPT: Compute cell indices once, reuse for all columns
cells <- cellFromXY(template, as.matrix(xy_all))

make_raster_fast <- function(values, template_rast, cells, name) {
  r <- template_rast
  r[] <- NA
  r[cells] <- values
  names(r) <- name
  r
}

.timer("Rasterise zones + uncertainty", {
  zone_rast   <- make_raster_fast(hard_zones,  template, cells, "zone")
  uncert_rast <- make_raster_fast(uncertainty, template, cells, "uncertainty")

  member_rasts <- lapply(seq_len(k_opt), function(j) {
    make_raster_fast(membership_mat[, j], template, cells, paste0("zone", j))
  })
  member_stack <- rast(member_rasts)
})

.timer("Write zone rasters", {
  writeRaster(zone_rast,    file.path(CFG$out_data, "zones_fuzzy.tif"),
              overwrite = TRUE)
  writeRaster(uncert_rast,  file.path(CFG$out_data, "zones_uncertainty.tif"),
              overwrite = TRUE)
  writeRaster(member_stack, file.path(CFG$out_data, "zones_membership.tif"),
              overwrite = TRUE)
})

# ── Zone + uncertainty maps ───────────────────────────────────────────────────
# Downsample rasters for plotting to avoid rendering millions of pixels
# in ggplot (which is very slow). Full-resolution rasters are saved above.

downsample_for_plot <- function(r, max_cells) {
  if (ncell(r) > max_cells) {
    fact <- ceiling(sqrt(ncell(r) / max_cells))
    aggregate(r, fact = fact, fun = "modal", na.rm = TRUE)
  } else {
    r
  }
}

downsample_continuous <- function(r, max_cells) {
  if (ncell(r) > max_cells) {
    fact <- ceiling(sqrt(ncell(r) / max_cells))
    aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
  } else {
    r
  }
}

zone_rast_plot   <- downsample_for_plot(zone_rast, CFG$plot_max_cells)
uncert_rast_plot <- downsample_continuous(uncert_rast, CFG$plot_max_cells)

zone_cols <- hcl.colors(k_opt, palette = "Set2")

p_zone <- ggplot() +
  geom_spatraster(data = as.factor(zone_rast_plot)) +
  scale_fill_manual(values    = setNames(zone_cols, as.character(1:k_opt)),
                    name      = "Zone",
                    na.value  = "grey90") +
  geom_sf(data = pts_sf_utm, colour = "white", size = 1.2, shape = 16) +
  labs(title    = "Management zones — fuzzy c-means",
       subtitle = paste0("k = ", k_opt, ", m = ", CFG$fuzziness,
                         "  |  n = ", nrow(soil_raw), " soil profiles")) +
  theme_bw(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

p_uncert <- ggplot() +
  geom_spatraster(data = uncert_rast_plot) +
  scale_fill_viridis_c(option   = "A",
                       name     = "Uncertainty\n(1 - membership)",
                       na.value = "grey90") +
  labs(title = "Zone assignment uncertainty") +
  theme_bw(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

p_zone_uncert <- (p_zone | p_uncert) + plot_annotation(tag_levels = "A")
ggsave(file.path(CFG$out_figures, "12_zone_maps.png"),
       p_zone_uncert, width = 13, height = 6, dpi = 300)

# ── Membership gradient maps ──────────────────────────────────────────────────
# Downsample before base-R plot as well
member_plot <- downsample_continuous(member_stack, CFG$plot_max_cells)

png(file.path(CFG$out_figures, "13_membership_maps.png"),
    width = 300 * k_opt, height = 600, res = 150)
plot(member_plot, col = viridis(100),
     main  = paste0("Membership — Zone ", 1:k_opt),
     axes  = FALSE, mar = c(1, 1, 2.5, 3), range = c(0, 1))
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 — SPATIALLY CONSTRAINED CLUSTERING (ClustGeo)
# ─────────────────────────────────────────────────────────────────────────────
# Use parallelDist for both D0 and D1 distance matrices (C++ threads).
# ─────────────────────────────────────────────────────────────────────────────

cat("\nRunning ClustGeo (k =", k_opt, ")...\n")
cat("Subsampling", CFG$sample_size_geo, "pixels for computation...\n")

set.seed(77)
idx_geo  <- sample(n_pixels, min(CFG$sample_size_geo, n_pixels))
sub_attr <- mat_scl[idx_geo, ]
sub_xy   <- as.matrix(xy_all[idx_geo, ]) / CFG$geo_unit

# parallelDist instead of base dist() — significant speedup for 3000 rows
.timer("ClustGeo distance matrices", {
  D0 <- parallelDist(sub_attr, method = "euclidean")
  D1 <- parallelDist(sub_xy,   method = "euclidean")
})

.timer("choicealpha", {
  alpha_choice <- choicealpha(
    D0          = D0,
    D1          = D1,
    range.alpha = CFG$alpha_grid,
    K           = k_opt,
    graph       = FALSE
  )
})

alpha_df <- tibble(
  alpha = CFG$alpha_grid,
  Q1    = alpha_choice$Qnorm[, 1],
  Q2    = alpha_choice$Qnorm[, 2]
)

p_alpha <- ggplot(alpha_df |> pivot_longer(-alpha), aes(alpha, value, colour = name)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(Q1 = "#3D7EBF", Q2 = "#E8813A"),
                      labels = c(Q1 = "Attribute homogeneity (Q1)",
                                 Q2 = "Spatial compactness (Q2)")) +
  scale_x_continuous(breaks = CFG$alpha_grid) +
  labs(title    = "ClustGeo: alpha selection",
       subtitle = "Choose alpha where Q2 >= 0.85 and Q1 is still acceptable",
       x = "alpha (mixing parameter)", y = "Normalised quality (Q)",
       colour = NULL) +
  geom_vline(xintercept = 0.2, linetype = "dashed", colour = "grey50") +
  theme_bw(base_size = 11)

ggsave(file.path(CFG$out_figures, "14_clustgeo_alpha.png"),
       p_alpha, width = 8, height = 4.5, dpi = 300)

alpha_opt <- alpha_df |>
  filter(Q2 >= 0.85, Q1 > 0.75) |>
  slice_max(Q1, n = 1) |>
  pull(alpha)

if (length(alpha_opt) == 0) {
  alpha_opt <- 0.2
  cat("Warning: no alpha meets both criteria; defaulting to alpha = 0.2\n")
}
cat("Selected alpha:", alpha_opt, "\n")
cat(">>> Inspect figures/14_clustgeo_alpha.png and adjust if needed.\n")

.timer("hclustgeo + cutree", {
  hclust_geo <- hclustgeo(D0, D1, alpha = alpha_opt)
  zones_geo  <- cutree(hclust_geo, k = k_opt)
})

df_geo_sub <- xy_all[idx_geo, ] |>
  as_tibble() |>
  mutate(zone_geo = zones_geo)

zone_geo_rast <- rast(
  data.frame(x = df_geo_sub$x, y = df_geo_sub$y, z = df_geo_sub$zone_geo),
  type = "xyz"
)
crs(zone_geo_rast) <- paste0("EPSG:", CFG$crs_code)
zone_geo_rast <- resample(zone_geo_rast, template, method = "near")

writeRaster(zone_geo_rast,
            file.path(CFG$out_data, "zones_clustgeo.tif"),
            overwrite = TRUE)

# OPT: Downsample for comparison plot
zone_geo_rast_plot <- downsample_for_plot(zone_geo_rast, CFG$plot_max_cells)

p_geo <- ggplot() +
  geom_spatraster(data = as.factor(zone_geo_rast_plot)) +
  scale_fill_manual(values   = setNames(zone_cols, as.character(1:k_opt)),
                    name     = "Zone",
                    na.value = "grey90") +
  labs(title    = "Management zones — ClustGeo",
       subtitle = paste("alpha =", alpha_opt)) +
  theme_bw(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

p_zones_comp <- (p_zone | p_geo) +
  plot_annotation(title = "Fuzzy c-means vs spatially constrained (ClustGeo)",
                  tag_levels = "A")
ggsave(file.path(CFG$out_figures, "15_zones_comparison.png"),
       p_zones_comp, width = 13, height = 6, dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6 — ZONE PROFILING WITH MEASURED SOIL PROPERTIES
# ─────────────────────────────────────────────────────────────────────────────

obs_zone <- terra::extract(zone_rast, vect(pts_sf_utm), ID = FALSE)[["zone"]]

soil_zones <- pts_sf_utm |>
  st_drop_geometry() |>
  select(all_of(target_vars), x, y) |>
  mutate(zone = paste0("Zone ", obs_zone)) |>
  drop_na(zone) |>
  filter(!is.na(obs_zone))

zone_summary <- soil_zones |>
  group_by(zone) |>
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

# ── Box plots ──────────────────────────────────────────────────────────────────
soil_long <- soil_zones |>
  pivot_longer(all_of(target_vars),
               names_to = "property", values_to = "value") |>
  mutate(property = factor(property, levels = target_vars))

ggplot(soil_long, aes(zone, value, fill = zone)) +
  geom_boxplot(outlier.shape = 21, linewidth = 0.6, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 18,
               size = 3, colour = "white") +
  facet_wrap(~property, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = setNames(zone_cols, paste0("Zone ", 1:k_opt))) +
  labs(title    = "Soil property distribution by management zone",
       subtitle = "Diamond = mean", x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

ggsave(file.path(CFG$out_figures, "16_zone_boxplots.png"),
       width = 11, height = 7, dpi = 300)

# ── Radar chart ────────────────────────────────────────────────────────────────
radar_df <- zone_summary |>
  select(zone, ends_with("_mean")) |>
  rename_with(~str_remove(.x, "_mean"), ends_with("_mean")) |>
  pivot_longer(-zone, names_to = "property", values_to = "value") |>
  group_by(property) |>
  mutate(norm = (value - min(value)) / (max(value) - min(value) + 1e-9)) |>
  ungroup()

ggplot(radar_df, aes(x = property, y = norm, colour = zone, group = zone)) +
  geom_polygon(aes(fill = zone), alpha = 0.15, linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = setNames(zone_cols, paste0("Zone ", 1:k_opt))) +
  scale_fill_manual(values   = setNames(zone_cols, paste0("Zone ", 1:k_opt))) +
  coord_polar() +
  labs(title = "Zone profile radar chart (normalised)",
       x = NULL, y = NULL, colour = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom")

ggsave(file.path(CFG$out_figures, "17_zone_radar.png"),
       width = 7, height = 7, dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 7 — FERTILITY PRESCRIPTIONS
# ─────────────────────────────────────────────────────────────────────────────

classify <- function(val, low, high, low_label, high_label,
                     mid_label = "Marginal") {
  case_when(
    val < low  ~ low_label,
    val > high ~ high_label,
    TRUE       ~ mid_label
  )
}

prescriptions <- zone_summary |>
  select(zone, n, ends_with("_mean")) |>
  rename_with(~str_remove(.x, "_mean")) |>
  mutate(
    SOC_status = if ("SOC" %in% names(pick(everything()))) classify(SOC,
                        CFG$thresh$SOC_low, CFG$thresh$SOC_high,
                        "Deficient — add organic amendments",
                        "Adequate") else NA_character_,

    pH_status  = if ("pH" %in% names(pick(everything()))) case_when(
      pH < CFG$thresh$pH_low  ~ "Acidic — lime application recommended",
      pH > CFG$thresh$pH_high ~ "Alkaline — gypsum or S may be needed",
      TRUE                    ~ "Optimal pH range"
    ) else NA_character_,

    CEC_status = if ("CEC" %in% names(pick(everything()))) classify(CEC,
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
  select(zone, n, any_of(target_vars),
         any_of(c("SOC_status", "pH_status", "CEC_status")), VR_priority)

cat("\n===  FERTILITY PRESCRIPTIONS  ===\n")
print(prescriptions)
write_csv(prescriptions,
          file.path(CFG$out_results, "fertility_prescriptions.csv"))

# ── Prescription summary bar chart ───────────────────────────────────────────
prescriptions |>
  mutate(zone = fct_reorder(zone, as.integer(str_extract(zone, "\\d+")))) |>
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
  ggplot(aes(zone, fill = severity)) +
  geom_bar(position = "stack") +
  facet_wrap(~indicator, nrow = 1) +
  scale_fill_manual(values = c("Action needed" = "#d73027",
                               "Monitor"       = "#fee090",
                               "Adequate"      = "#1a9850"),
                    name = NULL) +
  labs(title = "Soil fertility status by management zone",
       x = NULL, y = "Count of indicators") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")

ggsave(file.path(CFG$out_figures, "18_prescription_summary.png"),
       width = 10, height = 5, dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 8 — FINAL COMBINED FIGURE
# ─────────────────────────────────────────────────────────────────────────────

p_map_final <- ggplot() +
  geom_spatraster(data = as.factor(zone_rast_plot)) +
  scale_fill_manual(values   = setNames(zone_cols, as.character(1:k_opt)),
                    name     = "Zone",
                    na.value = "grey90") +
  geom_sf(data = pts_sf_utm, colour = "black", fill = "white",
          size = 2, shape = 21, stroke = 0.5) +
  labs(title    = "Management zones (fuzzy c-means)",
       subtitle = paste0("k = ", k_opt,
                         " | embedding PCs: 1–", n_pcs_cluster)) +
  theme_bw(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

box_var <- target_vars[1]
p_box <- soil_zones |>
  ggplot(aes(zone, .data[[box_var]], fill = zone)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 1.5) +
  scale_fill_manual(values = setNames(zone_cols, paste0("Zone ", 1:k_opt))) +
  labs(title    = paste(box_var, "by zone"),
       subtitle = "Dashed: deficiency / adequacy thresholds",
       x = NULL, y = box_var) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

(p_map_final | (p_box / p_alpha)) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(
    title = "Satellite Embedding Zones for Soil Fertility Assessment",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(file.path(CFG$out_figures, "19_final_figure.png"),
       width = 15, height = 8, dpi = 300)

cat("\n=== Script 04 complete ===\n")
cat("Outputs:\n")
cat("  Rasters  :", file.path(CFG$out_data,    "zones_*.tif"), "\n")
cat("  Tables   :", file.path(CFG$out_results, "*.csv"), "\n")
cat("  Figures  :", file.path(CFG$out_figures, "11-19_*.png"), "\n")

# =============================================================================
# END OF SCRIPT 04
# =============================================================================

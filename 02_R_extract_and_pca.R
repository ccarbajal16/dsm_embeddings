# =============================================================================
# Script 02: Extract embeddings at soil points + PCA on satellite embeddings
#
# Inputs:
#   data/soilgrids_250m/soilgrids_at_points.csv  — soil properties at points
#                                                   (cols: fid, longitude, latitude,
#                                                    BD, CEC, Fragm, Sand, Silt,
#                                                    Clay, N, OCD, pH, SOC)
#   data/embedding_raster_2024.tif               — 64-band GeoTIFF (A00…A63)
#
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
required <- c(
  "terra",       # raster / vector operations
  "sf",          # vector operations and CRS handling
  "tidyverse",   # data wrangling + ggplot2
  "factoextra",  # PCA visualisation helpers
  "umap",        # UMAP non-linear embedding
  "viridis",     # colour palettes
  "patchwork",   # compositing ggplot panels
  "tidyterra"    # ggplot2 + terra integration
)

new_pkgs <- required[!(required %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required, library, character.only = TRUE))

# ── 1. USER CONFIGURATION ─────────────────────────────────────────────────────

CFG <- list(
  # File paths
  soil_pts_csv   = "data/soilgrids_250m/soilgrids_at_points.csv",
  emb_raster     = "data/embedding_raster_2024.tif",

  # CRS of the embedding raster (EPSG:32717 — UTM Zone 17S, WGS84)
  crs_code       = 32717,

  # Coordinate columns in the soil CSV (WGS84 lon/lat)
  lon_col = "longitude",
  lat_col = "latitude",

  # Properties to predict (must match column names in soil CSV)
  target_vars = c("BD", "CEC", "Fragm", "Sand", "Silt", "Clay",
                  "N", "OCD", "pH", "SOC"),

  # PCA: number of components to retain (NULL = auto at 90% variance)
  n_pcs = NULL,

  # Output directories
  out_data    = "data/",
  out_figures = "figures/",
  out_results = "outputs/"
)

# Create output directories
walk(c(CFG$out_data, CFG$out_figures, CFG$out_results),
     ~if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# =============================================================================
# STEP 1 — EXTRACT EMBEDDING VALUES AT SOIL PROFILE LOCATIONS
# (logic from extract_values.R, integrated here)
# =============================================================================

cat("\n── Step 1: Extracting embedding values at soil points ────────────────────\n")

emb <- rast(CFG$emb_raster)
pts <- read.csv(CFG$soil_pts_csv)

cat("Soil points loaded:", nrow(pts), "rows\n")
cat("Raster CRS:", as.character(crs(emb, describe = TRUE)[["code"]]), "\n")
cat("Raster bands:", nlyr(emb), "\n")

# Project WGS84 points to the raster CRS
pts_sf <- st_as_sf(pts, coords = c(CFG$lon_col, CFG$lat_col), crs = 4326)
pts_sf <- st_transform(pts_sf, crs(emb))

# Extract 64-band embedding values
vals <- terra::extract(emb, vect(pts_sf))

# Combine original data with extracted values (drop the 'ID' column from extract)
result <- cbind(pts, vals[, -1])   # vals col 1 is the row ID from terra::extract

write.csv(result, "data/soilgrids_with_embeddings.csv", row.names = FALSE)
cat("Saved: data/soilgrids_with_embeddings.csv\n")

# =============================================================================
# MODULE 4 — PCA ON EMBEDDING BANDS
# =============================================================================

cat("\n── Module 4: PCA on 64-band satellite embeddings ─────────────────────────\n")

# ── 2. LOAD AND PREPARE SOIL DATA ─────────────────────────────────────────────

soil_raw <- read_csv("data/soilgrids_with_embeddings.csv", show_col_types = FALSE)

cat("\n── CSV column names ──────────────────────────────────────────────────────\n")
print(names(soil_raw))
cat("─────────────────────────────────────────────────────────────────────────\n\n")

# Auto-detect embedding columns (A00…A63 naming from GEE)
emb_cols <- grep("^A\\d{2}$",         names(soil_raw), value = TRUE)
if (!length(emb_cols))
  emb_cols <- grep("^embedding_\\d+$", names(soil_raw), value = TRUE)
if (!length(emb_cols))
  emb_cols <- grep("^B\\d+$",          names(soil_raw), value = TRUE)
if (!length(emb_cols))
  emb_cols <- grep("^emb_\\d+$",       names(soil_raw), value = TRUE)
if (!length(emb_cols))
  emb_cols <- grep("^first_\\d+$",     names(soil_raw), value = TRUE)

if (length(emb_cols) == 0) {
  stop(
    "\nNo embedding columns detected in the CSV.\n",
    "Look at the column names printed above and set emb_cols manually."
  )
}
cat("Embedding bands found:", length(emb_cols),
    "  (e.g.", emb_cols[1], "…", tail(emb_cols, 1), ")\n")

soil <- soil_raw

# Reproject points to UTM and add x/y columns for spatial operations
pts_soil_sf <- st_as_sf(soil, coords = c(CFG$lon_col, CFG$lat_col),
                         crs = 4326) |>
  st_transform(CFG$crs_code)

soil$x <- st_coordinates(pts_soil_sf)[, 1]
soil$y <- st_coordinates(pts_soil_sf)[, 2]

x_col <- "x"
y_col <- "y"

# Separate embedding matrix from metadata
emb_mat <- soil |> select(all_of(emb_cols)) |> as.matrix()
meta    <- soil |> select(all_of(c(x_col, y_col, CFG$target_vars)))

# Remove rows with any NA
complete_idx <- complete.cases(emb_mat) & complete.cases(meta)
emb_mat      <- emb_mat[complete_idx, ]
meta         <- meta[complete_idx, ]
cat("Complete profiles for analysis:", nrow(meta), "\n")

# ── 3. PCA ────────────────────────────────────────────────────────────────────

pca_obj  <- prcomp(emb_mat, center = TRUE, scale. = TRUE)
var_prop <- summary(pca_obj)$importance[2, ]   # proportion per PC
cum_var  <- summary(pca_obj)$importance[3, ]   # cumulative

# Auto-select n PCs explaining >= 90% variance (unless overridden)
n_pcs <- if (!is.null(CFG$n_pcs)) {
  CFG$n_pcs
} else {
  which(cum_var >= 0.90)[1]
}
cat("PCs retained (>=90% variance):", n_pcs, "\n")

# ── 3a. Scree plot ────────────────────────────────────────────────────────────
p_scree <- fviz_eig(pca_obj, ncp = 20, addlabels = TRUE,
                    barfill = "#3D7EBF", barcolor = NA) +
  geom_vline(xintercept = n_pcs + 0.5, linetype = "dashed",
             colour = "firebrick", linewidth = 0.8) +
  labs(title = "Scree plot — Satellite Embedding PCA",
       subtitle = sprintf("Dashed line: PC%d = %.1f%% cumulative variance",
                          n_pcs, cum_var[n_pcs] * 100)) +
  theme_bw(base_size = 11)

ggsave(file.path(CFG$out_figures, "01_scree_plot.png"),
       p_scree, width = 8, height = 4, dpi = 300)
cat("Saved: figures/01_scree_plot.png\n")

# ── 3b. PC scores tibble ──────────────────────────────────────────────────────
scores_df <- as_tibble(pca_obj$x[, 1:n_pcs]) |>
  bind_cols(meta)

# ── 3c. Biplots coloured by soil properties ───────────────────────────────────
make_biplot <- function(colour_var, palette = "D") {
  ggplot(scores_df, aes(PC1, PC2, colour = .data[[colour_var]])) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_colour_viridis_c(option = palette) +
    labs(
      title   = paste("PCA biplot —", colour_var),
      x       = sprintf("PC1 (%.1f%%)", var_prop[1] * 100),
      y       = sprintf("PC2 (%.1f%%)", var_prop[2] * 100),
      colour  = colour_var
    ) +
    theme_bw(base_size = 10)
}

tv       <- CFG$target_vars
n_tv     <- length(tv)
palettes <- rep(c("D", "C", "B", "A", "E"), length.out = n_tv)
bp_plots <- lapply(seq_len(min(n_tv, 4)),
                   function(i) make_biplot(tv[i], palettes[i]))
biplot_grid <- wrap_plots(bp_plots, ncol = 2) +
  plot_annotation(title = "PCA of 64-band Satellite Embedding",
                  tag_levels = "A")

ggsave(file.path(CFG$out_figures, "02_biplots.png"),
       biplot_grid, width = 12, height = 9, dpi = 300)
cat("Saved: figures/02_biplots.png\n")

# ── 3d. Spearman correlation heatmap: PC scores vs soil properties ─────────────
cor_df <- expand_grid(
  PC       = paste0("PC", 1:min(n_pcs, 15)),
  Property = CFG$target_vars
) |>
  rowwise() |>
  mutate(
    rho = cor(scores_df[[PC]], scores_df[[Property]],
              method = "spearman", use = "complete.obs"),
    p   = cor.test(scores_df[[PC]], scores_df[[Property]],
                   method = "spearman", exact = FALSE)$p.value
  ) |>
  ungroup() |>
  mutate(PC = factor(PC, levels = paste0("PC", 1:min(n_pcs, 15))))

p_cor <- ggplot(cor_df, aes(PC, Property, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf(ifelse(p < 0.05, "%.2f*", "%.2f"), rho)),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#d73027", mid = "white", high = "#1a9850",
    midpoint = 0, limits = c(-1, 1), name = "rho"
  ) +
  labs(title = "Spearman correlation: PC scores x soil properties",
       subtitle = "* p < 0.05", x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(CFG$out_figures, "03_pc_soil_correlation.png"),
       p_cor, width = 10, height = 5, dpi = 300)
cat("Saved: figures/03_pc_soil_correlation.png\n")

# ── 3e. UMAP visualisation ────────────────────────────────────────────────────
set.seed(42)
umap_res <- umap(emb_mat,
                 n_neighbors = 15, min_dist = 0.1,
                 n_components = 2, metric = "euclidean")

umap_df <- as_tibble(umap_res$layout) |>
  rename(UMAP1 = V1, UMAP2 = V2) |>
  bind_cols(meta)

make_umap <- function(colour_var, palette = "D") {
  ggplot(umap_df, aes(UMAP1, UMAP2, colour = .data[[colour_var]])) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_colour_viridis_c(option = palette) +
    labs(title = paste("UMAP —", colour_var)) +
    theme_bw(base_size = 10)
}

um_plots  <- lapply(seq_len(min(n_tv, 4)),
                    function(i) make_umap(tv[i], palettes[i]))
umap_grid <- wrap_plots(um_plots, ncol = 2) +
  plot_annotation(title = "UMAP of 64-band embedding space",
                  tag_levels = "A")

ggsave(file.path(CFG$out_figures, "04_umap.png"),
       umap_grid, width = 12, height = 9, dpi = 300)
cat("Saved: figures/04_umap.png\n")

# ── 3f. Project PCA onto the raster domain ───────────────────────────────────

cat("\nLoading embedding raster (may take a minute)...\n")
emb_rast <- rast(CFG$emb_raster)

cat("Projecting PCA onto raster grid...\n")
pca_rast <- predict(emb_rast, pca_obj, index = 1:n_pcs)
names(pca_rast) <- paste0("PC", 1:n_pcs)

writeRaster(pca_rast,
            file.path(CFG$out_data, "embedding_pca.tif"),
            overwrite = TRUE)
cat("Saved: data/embedding_pca.tif\n")

# Quick map of the first 4 PCs
png(file.path(CFG$out_figures, "05_pca_raster_maps.png"),
    width = 1600, height = 1200, res = 180)
plot(pca_rast[[1:min(4, n_pcs)]], col = viridis(100),
     main = paste0("PC", 1:min(4, n_pcs)), axes = FALSE, mar = c(1, 1, 2.5, 3))
dev.off()
cat("Saved: figures/05_pca_raster_maps.png\n")

# NOTE: pca_rast is NOT saved here — terra SpatRaster objects hold internal C++
# pointers that become invalid after save()/load() serialisation.
# Script 03 re-reads the raster directly from data/embedding_pca.tif.
save(CFG, scores_df, meta, emb_mat, emb_cols,
     x_col, y_col, n_pcs, pca_obj,
     file = "data/pca_workspace.RData")
cat("\nWorkspace saved: data/pca_workspace.RData\n")

cat("\n=== Script 02 complete ===\n")
cat("Outputs saved in:", CFG$out_data, CFG$out_figures, CFG$out_results, "\n")
cat("Next: run 03_R_spatial_prediction.R\n")

# =============================================================================
# END OF SCRIPT 02
# =============================================================================

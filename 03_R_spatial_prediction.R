# =============================================================================
# Script 03: Spatial Prediction of Soil Properties with Random Forest
#
# Inputs  (from Script 02):
#   data/pca_workspace.RData         — PCA objects, scores_df, meta, emb_mat,
#                                      n_pcs, pca_obj, pca_rast
#   data/soilgrids_with_embeddings.csv — soil data + embedding values at points
#
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
required <- c(
  "terra",          # raster operations
  "sf",             # vector operations
  "tidyverse",      # data wrangling + ggplot2
  "ranger",         # fast RF — used for cross-validation (Step 2)
  "caret",          # cross-validation helpers
  "randomForest",   # RF for spatial prediction (pure R, no Rcpp pointers)
  "viridis",        # colour palettes
  "patchwork",      # multi-panel plots
  "tidyterra"       # ggplot2 + terra integration
)

new_pkgs <- required[!(required %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required, library, character.only = TRUE))

# ── 1. USER CONFIGURATION ─────────────────────────────────────────────────────

CFG <- list(
  # Input files
  workspace      = "data/pca_workspace.RData",
  soil_csv       = "data/soilgrids_with_embeddings.csv",

  # CRS of the embedding raster (EPSG:32717 — UTM Zone 17S, WGS84)
  crs_code       = 32717,

  # Coordinate columns in the soil CSV (WGS84 lon/lat)
  lon_col        = "longitude",
  lat_col        = "latitude",

  # Soil properties to predict (must match column names in CSV)
  target_vars    = c("BD", "CEC", "Fragm", "Sand", "Silt", "Clay",
                     "N", "OCD", "pH", "SOC"),

  # Random Forest hyperparameters
  rf_trees       = 500,          # number of trees
  rf_mtry_frac   = 1/3,          # fraction of predictors tried at each split
  rf_min_node    = 5,            # minimum node size

  # Cross-validation
  cv_folds       = 5,            # number of CV folds
  cv_repeats     = 3,            # repeated k-fold repeats
  cv_seed        = 42,

  # Output directories
  out_data    = "data/",
  out_figures = "figures/",
  out_results = "outputs/"
)

walk(c(CFG$out_data, CFG$out_figures, CFG$out_results),
     ~if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# =============================================================================
# STEP 1 — LOAD PCA WORKSPACE FROM SCRIPT 02
# =============================================================================

cat("\n── Step 1: Loading PCA workspace ─────────────────────────────────────────\n")

# Preserve Script 03 CFG before loading workspace (load() would overwrite it)
CFG_03 <- CFG
load(CFG_03$workspace)
# Objects loaded from Script 02: CFG (02's), scores_df, meta, emb_mat,
#                 emb_cols, x_col, y_col, n_pcs, pca_obj
CFG <- CFG_03   # restore Script 03 configuration
rm(CFG_03)

# Re-read pca_rast from file — terra SpatRaster objects hold internal C++
# pointers that become invalid after save()/load(), so we never serialise them.
pca_rast <- rast(file.path(CFG$out_data, "embedding_pca.tif"))
names(pca_rast) <- paste0("PC", seq_len(nlyr(pca_rast)))

cat("PCs available:", n_pcs, "\n")
cat("Soil profiles:", nrow(meta), "\n")
cat("Target variables:", paste(CFG$target_vars, collapse = ", "), "\n")

# Reload soil CSV to get target columns aligned with complete_idx
soil_raw <- read_csv(CFG$soil_csv, show_col_types = FALSE)

# Re-derive complete_idx (same logic as Script 02)
emb_cols_local <- grep("^A\\d{2}$", names(soil_raw), value = TRUE)
emb_mat_raw    <- soil_raw |> select(all_of(emb_cols_local)) |> as.matrix()
meta_raw       <- soil_raw |> select(all_of(c(CFG$lon_col, CFG$lat_col,
                                               CFG$target_vars)))
complete_idx   <- complete.cases(emb_mat_raw) & complete.cases(meta_raw)
soil_complete  <- soil_raw[complete_idx, ]

cat("Complete profiles used for modelling:", sum(complete_idx), "\n")

# PC scores at point locations (from prcomp object)
pc_scores <- as_tibble(pca_obj$x[, 1:n_pcs])
names(pc_scores) <- paste0("PC", 1:n_pcs)

# Bind target variables
model_df <- bind_cols(
  pc_scores,
  soil_complete |> select(all_of(CFG$target_vars))
)

# =============================================================================
# STEP 2 — CROSS-VALIDATED RANDOM FOREST
# =============================================================================

cat("\n── Step 2: Fitting RF models (", CFG$cv_folds, "×", CFG$cv_repeats,
    "CV) ──────────────────────────────\n", sep = "")

set.seed(CFG$cv_seed)
cv_ctrl <- trainControl(
  method      = "repeatedcv",
  number      = CFG$cv_folds,
  repeats     = CFG$cv_repeats,
  savePredictions = "final",
  allowParallel   = FALSE
)

pc_predictors <- paste0("PC", 1:n_pcs)
rf_models     <- list()
cv_results    <- list()

for (tvar in CFG$target_vars) {

  cat("  Fitting RF for:", tvar, "... ")

  # Drop rows where the target is NA
  df_tvar <- model_df |>
    select(all_of(c(pc_predictors, tvar))) |>
    drop_na()

  if (nrow(df_tvar) < 20) {
    cat("skipped (n =", nrow(df_tvar), "< 20)\n")
    next
  }

  mtry_val <- max(1, round(n_pcs * CFG$rf_mtry_frac))

  suppressWarnings(
    fit <- train(
      x          = df_tvar[, pc_predictors],
      y          = df_tvar[[tvar]],
      method     = "ranger",
      trControl  = cv_ctrl,
      tuneGrid   = data.frame(
        mtry              = mtry_val,
        splitrule         = "variance",
        min.node.size     = CFG$rf_min_node
      ),
      num.trees  = CFG$rf_trees,
      importance = "impurity"
    )
  )

  rf_models[[tvar]] <- fit

  # CV performance
  pred_obs <- fit$pred |>
    group_by(rowIndex) |>
    summarise(obs  = mean(obs),
              pred = mean(pred),
              .groups = "drop")

  rmse_val <- sqrt(mean((pred_obs$obs - pred_obs$pred)^2))
  mae_val  <- mean(abs(pred_obs$obs - pred_obs$pred))
  r2_val   <- cor(pred_obs$obs, pred_obs$pred, use = "complete.obs")^2

  cv_results[[tvar]] <- tibble(
    variable = tvar,
    n        = nrow(df_tvar),
    RMSE     = round(rmse_val, 4),
    MAE      = round(mae_val,  4),
    R2       = round(r2_val,   4)
  )

  cat("RMSE =", round(rmse_val, 3), " R² =", round(r2_val, 3), "\n")
}

cv_table <- bind_rows(cv_results)
write_csv(cv_table, file.path(CFG$out_results, "rf_cv_results.csv"))
cat("\nCV performance saved: outputs/rf_cv_results.csv\n")
print(cv_table)

# =============================================================================
# STEP 3 — VARIABLE IMPORTANCE
# =============================================================================

cat("\n── Step 3: Variable importance ────────────────────────────────────────────\n")

imp_list <- map(names(rf_models), function(tvar) {
  imp <- varImp(rf_models[[tvar]])$importance
  tibble(
    variable  = tvar,
    PC        = rownames(imp),
    Importance = imp[[1]]
  )
})

imp_df <- bind_rows(imp_list) |>
  mutate(PC = factor(PC, levels = paste0("PC", 1:n_pcs)))

write_csv(imp_df, file.path(CFG$out_results, "rf_importance.csv"))

# Importance heat-map
p_imp <- ggplot(imp_df, aes(PC, variable, fill = Importance)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_viridis_c(option = "C", name = "Importance\n(impurity)") +
  labs(title = "Random Forest variable importance",
       subtitle = paste0("n_pcs = ", n_pcs,
                         " | trees = ", CFG$rf_trees),
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(CFG$out_figures, "07_rf_importance.png"),
       p_imp, width = 10, height = 5, dpi = 300)
cat("Saved: figures/07_rf_importance.png\n")

# =============================================================================
# STEP 4 — OBSERVED vs PREDICTED SCATTER PLOTS
# =============================================================================

cat("\n── Step 4: Observed vs predicted plots ────────────────────────────────────\n")

obs_pred_list <- map(names(rf_models), function(tvar) {
  rf_models[[tvar]]$pred |>
    group_by(rowIndex) |>
    summarise(obs = mean(obs), pred = mean(pred), .groups = "drop") |>
    mutate(variable = tvar)
})
obs_pred_df <- bind_rows(obs_pred_list)

p_obs_pred <- ggplot(obs_pred_df, aes(obs, pred)) +
  geom_point(alpha = 0.5, size = 1.8, colour = "#3D7EBF") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "firebrick", linewidth = 0.8) +
  facet_wrap(~variable, scales = "free", ncol = 4) +
  labs(title  = "Observed vs predicted — cross-validated RF",
       subtitle = paste0(CFG$cv_folds, "-fold × ", CFG$cv_repeats, " repeats"),
       x = "Observed", y = "Predicted") +
  theme_bw(base_size = 10)

n_tvars <- length(names(rf_models))
fig_h <- ceiling(n_tvars / 4) * 3.5 + 1.5

ggsave(file.path(CFG$out_figures, "06_rf_obs_pred.png"),
       p_obs_pred, width = 14, height = fig_h, dpi = 300)
cat("Saved: figures/06_rf_obs_pred.png\n")

# =============================================================================
# STEP 5 — PROJECT PREDICTIONS ONTO THE RASTER DOMAIN
# =============================================================================

cat("\n── Step 5: Projecting RF predictions onto raster ─────────────────────────\n")
cat("PCA raster bands:", n_pcs, "\n")

# NOTE: ranger uses Rcpp C++ external pointers that can become invalid
# ("external pointer is not valid" / "NULL value passed as symbol address")
# even in freshly-fit models when large data operations follow.
# Solution: use randomForest for spatial prediction — it is a pure R object
# with no external pointers. ranger is kept only for the faster CV in Step 2.

# Extract all raster pixel values once (keep NA rows so indices align with raster)
pca_df     <- as.data.frame(pca_rast, xy = FALSE, na.rm = FALSE)
names(pca_df) <- paste0("PC", seq_len(ncol(pca_df)))
valid_rows <- complete.cases(pca_df)
pca_valid  <- pca_df[valid_rows, , drop = FALSE]   # subset used for prediction
cat("Raster pixels to predict:", sum(valid_rows), "of", nrow(pca_df), "\n")

template_r <- pca_rast[[1]]   # template for rebuilding result rasters

pred_rasts <- list()

for (tvar in names(rf_models)) {
  cat("  Predicting:", tvar, "... ")

  df_tvar <- as.data.frame(
    model_df |> select(all_of(c(pc_predictors, tvar))) |> drop_na()
  )

  # Fit final model with randomForest (pure R, no Rcpp, reliable serialisation)
  final_rf <- randomForest::randomForest(
    x     = df_tvar[, pc_predictors, drop = FALSE],
    y     = df_tvar[[tvar]],
    ntree = CFG$rf_trees,
    mtry  = max(1L, round(n_pcs * CFG$rf_mtry_frac))
  )

  pred_vals             <- rep(NA_real_, nrow(pca_df))
  pred_vals[valid_rows] <- predict(final_rf, newdata = pca_valid)

  pred_r         <- template_r
  values(pred_r) <- pred_vals
  names(pred_r)  <- tvar
  pred_rasts[[tvar]] <- pred_r
  cat("done\n")
}

pred_stack <- rast(pred_rasts)   # terra::rast() stacks a list of SpatRasters

writeRaster(pred_stack,
            file.path(CFG$out_data, "soil_predictions.tif"),
            overwrite = TRUE)
cat("Saved: data/soil_predictions.tif\n")

# =============================================================================
# STEP 6 — PREDICTION MAPS
# =============================================================================

cat("\n── Step 6: Prediction maps ────────────────────────────────────────────────\n")

n_pred <- nlyr(pred_stack)
n_cols <- min(4, n_pred)
n_rows <- ceiling(n_pred / n_cols)

png(file.path(CFG$out_figures, "08_soil_prediction_maps.png"),
    width  = n_cols * 400,
    height = n_rows * 350 + 100,
    res    = 150)
plot(pred_stack,
     col  = viridis(100),
     main = names(pred_stack),
     axes = FALSE,
     mar  = c(1, 1, 2.5, 3))
dev.off()
cat("Saved: figures/08_soil_prediction_maps.png\n")

# =============================================================================
# SAVE WORKSPACE FOR SCRIPT 04
# =============================================================================

save(CFG, pred_stack, cv_table, rf_models, imp_df,
     file = "data/rf_workspace.RData")
cat("\nWorkspace saved: data/rf_workspace.RData\n")

cat("\n=== Script 03 complete ===\n")
cat("Outputs:\n")
cat("  Rasters  :", file.path(CFG$out_data,    "soil_predictions.tif"), "\n")
cat("  Tables   :", file.path(CFG$out_results, "rf_cv_results.csv"),
    "&", file.path(CFG$out_results, "rf_importance.csv"), "\n")
cat("  Figures  :", file.path(CFG$out_figures, "06–08_*.png"), "\n")

# =============================================================================
# END OF SCRIPT 03
# =============================================================================

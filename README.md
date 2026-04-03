# DSM with Satellite Embeddings

Digital Soil Mapping (DSM) of 10 soil properties using Google's 64-band satellite embedding dataset (`GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL`) as covariates, reduced via PCA, modelled with Random Forest, and delineated into management zones via fuzzy c-means and spatially constrained clustering.

## Overview

This workflow downloads a pre-computed satellite embedding raster from Google Earth Engine, extracts the embedding values at soil profile locations, reduces the 64 bands to a compact set of principal components, fits cross-validated Random Forest models to predict soil properties across the landscape, and delineates management zones using the embedding PCs for precision agriculture applications.

**Study area:** Trujillo region, Peru (UTM Zone 17S, EPSG:32717)
**Target year:** 2024
**Soil depth:** 15–30 cm (SoilGrids 250 m)

## Workflow

```
01_download_raster.js   →   02_R_extract_and_pca.R   →   03_R_spatial_prediction.R   →   04_R_management_zones.R
      (GEE)                          (R)                           (R)                               (R)
```

### Script 01 — Download embedding raster (Google Earth Engine)

Loads the `GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL` image collection, filters to the study area and target year, mosaics the tiles, clips to the boundary, and exports a 64-band Cloud-Optimised GeoTIFF to Google Drive.

**Configure before running:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `boundaryAsset` | `projects/ee-cmcarbajal/assets/border_trujillo` | GEE asset path for study area boundary |
| `year` | `2024` | Target year (2017–2025) |
| `crs` | `EPSG:32717` | Output CRS |
| `scale` | `10` | Export resolution in metres |
| `driveFolder` | `DSM_embeddings` | Google Drive output folder |

**Output:** `embedding_raster_<year>.tif` (64 bands: A00–A63) in Google Drive.

---

### Script 02 — Extract embeddings & PCA (R)

1. Extracts 64-band embedding values at each soil profile location
2. Runs PCA on the embedding matrix (auto-selects PCs explaining ≥ 90% cumulative variance)
3. Generates exploratory figures (scree plot, biplots, UMAP, PC–soil correlation heatmap)
4. Projects the PCA model onto the full raster grid and saves `embedding_pca.tif`

**Key output:** 15 PCs retained (90.2% cumulative variance)

---

### Script 03 — Spatial prediction with Random Forest (R)

1. Loads PCA workspace from Script 02
2. Fits RF models (via `ranger`) for each soil property using repeated k-fold CV (5-fold × 3 repeats)
3. Extracts variable importance and plots observed vs predicted
4. Re-fits final models with `randomForest` (pure R, avoids `ranger`/Rcpp serialisation issues) and predicts across the raster domain
5. Saves a 10-band prediction raster `soil_predictions.tif`

**RF hyperparameters:**
| Parameter | Value |
|-----------|-------|
| Trees | 500 |
| mtry fraction | 1/3 of PCs |
| Min node size | 5 |
| CV folds × repeats | 5 × 3 |

---

### Script 04 — Management Zone Delineation (R)

Uses the embedding PCs from Script 02 to delineate precision-agriculture management zones across the landscape.

1. **Zone number selection** — evaluates k = 2–10 using elbow (WSS) and average silhouette width via mini-batch k-means (`ClusterR`); auto-selects the silhouette peak or accepts a user-supplied `k_final`
2. **Fuzzy c-means (FCM)** — fits FCM (via `ppclust`) on a 30 000-pixel subsample, assigns full-grid memberships analytically, and writes three rasters: hard zones, uncertainty (1 − max membership), and per-zone membership gradients
3. **Spatially constrained clustering (ClustGeo)** — builds attribute (`D0`) and geographic (`D1`) distance matrices with `parallelDist`, sweeps `alpha` to balance attribute homogeneity vs spatial compactness, and performs Ward hierarchical clustering on the blend
4. **Zone profiling** — extracts measured soil properties at profile locations per zone, computes summary statistics, and produces box plots and a radar chart
5. **Fertility prescriptions** — classifies each zone for SOC, pH, and CEC against agronomic thresholds and assigns a variable-rate application priority (HIGH / MEDIUM / LOW)

**Key configuration (`CFG`):**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_pcs_cluster` | 10 | Number of PCs used for clustering |
| `k_min` / `k_max` | 2 / 10 | Search range for optimal k |
| `k_final` | `NULL` | Override auto-selection if needed |
| `fuzziness` | 2.0 | FCM fuzziness exponent (m) |
| `sample_size_fcm` | 30 000 | Pixels used to fit FCM |
| `alpha_grid` | 0–1 by 0.1 | ClustGeo mixing parameter sweep |

**Outputs:**

| File | Description |
|------|-------------|
| `data/zones_fuzzy.tif` | Hard management zone raster |
| `data/zones_uncertainty.tif` | Assignment uncertainty raster |
| `data/zones_membership.tif` | Per-zone fuzzy membership stack |
| `data/zones_clustgeo.tif` | Spatially constrained zone raster |
| `outputs/zone_profiles.csv` | Mean / SD / median per zone × soil property |
| `outputs/fertility_prescriptions.csv` | Fertility status and VR priority per zone |

---

## Predicted Soil Properties

| Variable | Description |
|----------|-------------|
| BD | Bulk density (g/cm³) |
| CEC | Cation exchange capacity (cmol/kg) |
| Fragm | Coarse fragments (%) |
| Sand | Sand content (%) |
| Silt | Silt content (%) |
| Clay | Clay content (%) |
| N | Total nitrogen (g/kg) |
| OCD | Organic carbon density (kg/m³) |
| pH | Soil pH (H₂O) |
| SOC | Soil organic carbon (g/kg) |

## Repository Structure

```
dsm_embeddings/
├── 01_download_raster.js          # GEE script — export embedding raster
├── 02_R_extract_and_pca.R         # R — extract values, PCA, UMAP
├── 03_R_spatial_prediction.R      # R — RF modelling, spatial prediction
├── 04_R_management_zones.R        # R — management zone delineation & prescriptions
│
├── data/
│   ├── agro_geo.gpkg              # Study area boundary (GeoPackage)
│   ├── soilgrids_250m/
│   │   ├── soilgrids_at_points.csv        # Soil properties at profile locations
│   │   ├── soilgrids_stack_utm.tif        # SoilGrids raster stack (UTM)
│   │   └── soilgrids_<property>_15-30cm.tif  # Individual SoilGrids bands
│   ├── soilgrids_with_embeddings.csv      # Soil data + extracted embedding values
│   ├── embedding_pca.tif          # PCA raster (n_pcs bands) [gitignored]
│   ├── embedding_raster_2024.tif  # Raw 64-band embedding raster [gitignored]
│   ├── pca_workspace.RData        # PCA objects for Script 03
│   ├── rf_workspace.RData         # RF models and results
│   ├── zones_fuzzy.tif            # Hard management zone raster
│   ├── zones_uncertainty.tif      # FCM assignment uncertainty raster
│   ├── zones_membership.tif       # Per-zone fuzzy membership stack
│   └── zones_clustgeo.tif         # Spatially constrained zone raster
│
├── outputs/
│   ├── rf_cv_results.csv          # Cross-validated RMSE, MAE, R² per property
│   ├── rf_importance.csv          # Variable importance per property × PC
│   ├── zone_profiles.csv          # Mean / SD / median per zone × soil property
│   └── fertility_prescriptions.csv # Fertility status and VR priority per zone
│
└── figures/
    ├── 01_scree_plot.png           # PCA variance explained
    ├── 02_biplots.png              # PC1 vs PC2 coloured by soil properties
    ├── 03_pc_soil_correlation.png  # Spearman ρ heatmap: PCs × soil properties
    ├── 04_umap.png                 # UMAP of 64-band embedding space
    ├── 05_pca_raster_maps.png      # Spatial maps of first 4 PCs
    ├── 06_rf_obs_pred.png          # Observed vs predicted (CV) per property
    ├── 07_rf_importance.png        # RF variable importance heatmap
    ├── 08_soil_prediction_maps.png # Final soil prediction maps
    ├── 11_zone_selection.png       # Elbow and silhouette curves
    ├── 12_zone_maps.png            # FCM zone map + uncertainty
    ├── 13_membership_maps.png      # Per-zone fuzzy membership gradients
    ├── 14_clustgeo_alpha.png       # ClustGeo alpha selection
    ├── 15_zones_comparison.png     # FCM vs ClustGeo comparison
    ├── 16_zone_boxplots.png        # Soil properties by zone (box plots)
    ├── 17_zone_radar.png           # Zone profile radar chart
    ├── 18_prescription_summary.png # Fertility status by zone
    └── 19_final_figure.png         # Combined zone map + box plot + alpha plot
```

## R Dependencies

All scripts auto-install any missing packages on first run.

| Package | Used in | Purpose |
|---------|---------|---------|
| `terra` | 02–04 | Raster and vector operations |
| `sf` | 02–04 | Vector geometry and CRS handling |
| `tidyverse` | 02–04 | Data wrangling and ggplot2 |
| `factoextra` | 02, 04 | PCA / cluster visualisation helpers |
| `umap` | 02 | UMAP non-linear dimensionality reduction |
| `ranger` | 03 | Fast Random Forest (cross-validation) |
| `randomForest` | 03 | RF for spatial prediction (serialisation-safe) |
| `caret` | 03 | Cross-validation helpers |
| `viridis` | 02–04 | Colour palettes |
| `patchwork` | 02–04 | Multi-panel plot composition |
| `tidyterra` | 02–04 | ggplot2 + terra integration |
| `ppclust` | 04 | Fuzzy c-means clustering |
| `cluster` | 04 | Silhouette index |
| `ClustGeo` | 04 | Spatially constrained hierarchical clustering |
| `ClusterR` | 04 | Mini-batch k-means (C++ backend) |
| `matrixStats` | 04 | Fast row-wise operations |
| `parallelDist` | 04 | Multi-threaded distance matrices (C++) |

## Outputs

### Figures

| Figure | Description |
|--------|-------------|
| ![Scree plot](figures/01_scree_plot.png) | 15 PCs explain 90.2% of embedding variance |
| ![Biplots](figures/02_biplots.png) | PC space coloured by BD, CEC, Fragm, Sand |
| ![PCA maps](figures/05_pca_raster_maps.png) | Spatial structure of first 4 PCs |
| ![Predictions](figures/08_soil_prediction_maps.png) | Predicted maps for all 10 soil properties |
| ![Management zones](figures/19_final_figure.png) | Management zones, soil box plot, and ClustGeo alpha selection |

### Key files

- `data/soil_predictions.tif` — 10-band GeoTIFF with one layer per soil property
- `outputs/rf_cv_results.csv` — RMSE, MAE, R² for each property from repeated CV
- `data/zones_fuzzy.tif` — Hard management zone raster (FCM)
- `outputs/zone_profiles.csv` — Soil property statistics per zone
- `outputs/fertility_prescriptions.csv` — Variable-rate application priorities per zone

## Quick Start

1. **Run Script 01** in the [Google Earth Engine Code Editor](https://code.earthengine.google.com/), submit the export task from the Tasks panel, and download the resulting GeoTIFF to `data/`.

2. **Run Script 02** in R from the project root:
   ```r
   source("02_R_extract_and_pca.R")
   ```

3. **Run Script 03** in R:
   ```r
   source("03_R_spatial_prediction.R")
   ```

4. **Run Script 04** in R:
   ```r
   source("04_R_management_zones.R")
   ```
   > Inspect `figures/11_zone_selection.png` to validate the chosen `k`. Override `CFG$k_final` if needed and re-run. Similarly, check `figures/14_clustgeo_alpha.png` to adjust the ClustGeo `alpha` parameter.

> All paths are relative to the project root. Set your working directory to the repo root before running R scripts.

## Notes

- Large raster files (`embedding_raster_2024.tif`, `embedding_pca.tif`) are excluded from version control via `.gitignore` (> 100 MB).
- Script 03 uses `randomForest` (not `ranger`) for pixel-level prediction to avoid Rcpp external pointer issues that arise when predicting over large rasters after serialisation.
- The embedding dataset covers years 2017–2025. Change `CONFIG.year` in Script 01 to use a different year.

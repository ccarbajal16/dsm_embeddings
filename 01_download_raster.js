// =============================================================================
// Script 01: Download Satellite Embedding Raster with GEE
//
// Dataset: GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL
//

// -- USER CONFIGURATION -------------------------------------------------------

var CONFIG = {
  // GEE asset path for the study area boundary
  boundaryAsset : 'projects/ee-cmcarbajal/assets/border_trujillo',

  // Target year for the annual embedding composite (valid range: 2017-2025)
  year : 2024,

  // Output CRS (UTM zone for your study area; change EPSG as needed)
  crs : 'EPSG:32717',

  // Export scale in metres (10 = native resolution;
  // increase to 30 or 100 to reduce file size for large areas)
  scale : 10,

  // Google Drive folder for the export
  driveFolder : 'DSM_embeddings',

  // Maximum number of pixels for image export (increase if export fails)
  maxPixels : 1e10
};

// -- CONFIG VALIDATION --------------------------------------------------------
if (CONFIG.year < 2017 || CONFIG.year > 2025) {
  throw new Error(
    'CONFIG.year must be between 2017 and 2025. Got: ' + CONFIG.year
  );
}

// =============================================================================
// STEP 1 - LOAD STUDY BOUNDARY
// =============================================================================

var boundary = ee.FeatureCollection(CONFIG.boundaryAsset);
var region   = boundary.geometry();

Map.centerObject(boundary);
Map.addLayer(boundary, {color: '0000FF', fillColor: '00000000'}, 'Study boundary');

// =============================================================================
// STEP 2 - LOAD THE SATELLITE EMBEDDING COLLECTION
// =============================================================================

var startDate = ee.Date.fromYMD(CONFIG.year,     1, 1);
var endDate   = ee.Date.fromYMD(CONFIG.year + 1, 1, 1);

var embeddingCollection = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
  .filterDate(startDate, endDate)
  .filterBounds(region);

// Guard: warn if no tile found for the requested year / region
embeddingCollection.size().evaluate(function(n) {
  if (n === 0) {
    print('WARNING: No embedding tile found for year ' + CONFIG.year +
          ' in the study area. Try a different year or check the boundary asset.');
  } else {
    print('Embedding tiles found for ' + CONFIG.year + ': ' + n);
  }
});

// Mosaic all tiles that cover the study area into one image
var embeddingImage = embeddingCollection.mosaic();

// Verify band count (should be 64: A00 ... A63)
var bandNames = embeddingImage.bandNames();
print('Embedding band count:', bandNames.length());
print('Band names (first 5):', bandNames.slice(0, 5));
print('Projection:', embeddingImage.select([0]).projection());

// =============================================================================
// STEP 3 - CLIP RASTER TO STUDY BOUNDARY
// =============================================================================

// updateMask is more memory-efficient than clip() for irregular shapes;
// pixels outside the boundary become no-data rather than clipped to bbox.
var embeddingClipped = embeddingImage.updateMask(
  ee.Image.constant(1).clip(region)
);

// Quick visualisation: RGB composite from first 3 bands
Map.addLayer(
  embeddingClipped.select(bandNames.slice(0, 3)),
  {min: -3, max: 3, gamma: 1.2},
  'Embedding bands 0-2 (RGB)'
);

// =============================================================================
// STEP 4 - EXPORT EMBEDDING RASTER (GeoTIFF)
// =============================================================================

var rasterName = 'embedding_raster_' + CONFIG.year;

Export.image.toDrive({
  image          : embeddingClipped,
  description    : rasterName,
  folder         : CONFIG.driveFolder,
  fileNamePrefix : rasterName,
  region         : region,
  scale          : CONFIG.scale,
  crs            : CONFIG.crs,
  maxPixels      : CONFIG.maxPixels,
  fileFormat     : 'GeoTIFF',
  formatOptions  : { cloudOptimized: true }
});

print('>>> Export task ready: ' + rasterName + '.tif');
print('Go to the Tasks panel (top right) and click Run to start the export.');
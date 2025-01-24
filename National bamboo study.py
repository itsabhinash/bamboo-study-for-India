# Bamboo study for india

# Geemap

import ee
import geemap
import json
import os
import requests
from geemap import geojson_to_ee, ee_to_geojson
from ipyleaflet import GeoJSON, Marker, MarkerCluster
import matplotlib.pyplot as plt
import pandas as pd
import numpy as num

ee.Initialize()

# Load the India Boundary (Study area)
India_boundary = ee.FeatureCollection('path/India_boundary')


# Tree presence
# ESA worldcover 10m
esa_lc = ee.ImageCollection('ESA/WorldCover/v200').first()

# Define the land cover classes to be masked.
classes_to_mask = [20,30,40,50,60,70,80,90,95,100]

# Create a mask where the undesired land cover classes are assigned a value of 0.
mask = esa_lc.neq(classes_to_mask[0])
for i in range(1, len(classes_to_mask)):
    mask = mask.And(esa_lc.neq(classes_to_mask[i]))

# Apply the mask to the land cover dataset to get tree presence areas.
masked_lc = esa_lc.mask(mask).clip(India_boundary)


# Processing for spectral data extraction
# Considered time periods
    #  Summer median = '2023-03-01 to '2023-05-31'
    #  Winter median = '2023-10-01 to '2023-12-31'
# Considered bands
    #
def maskCloudAndShadows(image):
    cloudProb = image.select('MSK_CLDPRB')
    snowProb = image.select('MSK_SNWPRB')
    cloud = cloudProb.lt(5)
    snow = snowProb.lt(5)
    scl = image.select('SCL')
    shadow = scl.eq(3)  # 3 = cloud shadow
    cirrus = scl.eq(10)  # 10 = cirrus
    # Cloud probability less than 5% or cloud shadow classification
    mask = (cloud.And(snow)).And(cirrus.neq(1)).And(shadow.neq(1))
    return image.updateMask(mask)

startDate = '2023-03-01'
endDate = '2023-05-31'

collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterDate(startDate, endDate).map(maskCloudAndShadows).filterBounds(India_boundary)
collection_median = collection.select('B1').median().updateMask(masked_lc).divide(10000).rename('collection_median')

# Exporting band data to drive
geemap.ee_export_image_to_drive(
   collection_median.clip(India_boundary) , 
   description="B1_India_summer", 
   folder="gee",
   scale=10 , 
   region=India_boundary.geometry().dissolve() , 
   maxPixels=1e13
)


# Processing for backscatter data extraction
# Considered time periods
    #  Summer median = '2023-03-01 to '2023-05-31'
    #  Winter median = '2023-10-01 to '2023-12-31'
# Considered bands
    # VV,VH

# Import the Sentinel-1 ImageCollection
sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')

# Function to filter Sentinel-1 data based on orbit and polarization pair
def get_filtered_images(orbit_pass, polarization_pair):
    return sentinel1.filterBounds(India_boundary) \
                    .filterDate('2023-03-01', '2023-05-31') \
                    .filter(ee.Filter.eq('orbitProperties_pass', orbit_pass)) \
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', polarization_pair[0])) \
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', polarization_pair[1])) \
                    .filter(ee.Filter.eq('instrumentMode', 'IW'))

# Filter for VV + VH polarizations
vv_vh_desc = get_filtered_images('DESCENDING', ['VV', 'VH'])

# Create median composites for each combination
combined = vv_vh_desc.select('VV', 'VH').median().clip(India_boundary).rename(['VV_desc', 'VH_desc']).updateMask(masked_lc)

# Exporting band data to drive
geemap.ee_export_image_to_drive(
    combined.clip(India_boundary),
    description="backscatter_summer",
    folder="gee",
    scale=10,
    region=India_boundary.geometry(),  
    maxPixels=1e13  
)


# Processing for soil data extraction
# loading data
bdod = ee.Image("projects/soilgrids-isric/bdod_mean") #Bulk density of the fine earth fraction
cec = ee.Image("projects/soilgrids-isric/cec_mean") #Cation Exchange Capacity of the soil	
cfvo =  ee.Image("projects/soilgrids-isric/cfvo_mean") #Volumetric fraction of coarse fragments (> 2 mm)
clay = ee.Image("projects/soilgrids-isric/clay_mean") #Proportion of clay particles (< 0.002 mm) in the fine earth fraction
nitrogen = ee.Image("projects/soilgrids-isric/nitrogen_mean") #Total nitrogen (N)
phh2o = ee.Image("projects/soilgrids-isric/phh2o_mean") #Soil pH
sand = ee.Image("projects/soilgrids-isric/sand_mean") #Proportion of sand particles (> 0.05 mm) in the fine earth fraction
silt = ee.Image("projects/soilgrids-isric/silt_mean") #Proportion of silt particles (≥ 0.002 mm and ≤ 0.05 mm) in the fine earth fraction
soc = ee.Image("projects/soilgrids-isric/soc_mean") #Soil organic carbon content in the fine earth fraction
ocd = ee.Image("projects/soilgrids-isric/ocd_mean") #Organic carbon density
ocs = ee.Image("projects/soilgrids-isric/ocs_mean") #Organic carbon stocks

# Preprocessing of each paramter
bdod_f = bdod.expression(
    '(0.005 * (0.5*((   (5*(ZERO+ONE)) + (10*(ONE+TWO)) + (15*(TWO + THREE))  +  (30*(THREE + FOUR))  +  (40*(FOUR + FIVE))  +   (100*(FIVE))    ) )  )) ', {
      'ZERO': bdod.select('bdod_0-5cm_mean'),
      'ONE': bdod.select('bdod_5-15cm_mean'),
      'TWO': bdod.select('bdod_15-30cm_mean'),
      'THREE': bdod.select('bdod_30-60cm_mean'),
      'FOUR': bdod.select('bdod_60-100cm_mean'),
      'FIVE': bdod.select('bdod_100-200cm_mean')
}).clip(India_boundary).divide(100).rename('bdod_f')

cec_f = cec.expression(
    '(0.005 * (0.5*((   (5*(ZERO+ONE)) + (10*(ONE+TWO)) + (15*(TWO + THREE))  +  (30*(THREE + FOUR))  +  (40*(FOUR + FIVE))  +   (100*(FIVE))    ) )  )) ', {
      'ZERO': cec.select('cec_0-5cm_mean'),
      'ONE': cec.select('cec_5-15cm_mean'),
      'TWO': cec.select('cec_15-30cm_mean'),
      'THREE': cec.select('cec_30-60cm_mean'),
      'FOUR': cec.select('cec_60-100cm_mean'),
      'FIVE': cec.select('cec_100-200cm_mean')
}).clip(India_boundary).divide(10).rename('cec_f')

cfvo_f = cfvo.expression(
    '(0.005 * (0.5*((   (5*(ZERO+ONE)) + (10*(ONE+TWO)) + (15*(TWO + THREE))  +  (30*(THREE + FOUR))  +  (40*(FOUR + FIVE))  +   (100*(FIVE))    ) )  )) ', {
      'ZERO': cfvo.select('cfvo_0-5cm_mean'),
      'ONE': cfvo.select('cfvo_5-15cm_mean'),
      'TWO': cfvo.select('cfvo_15-30cm_mean'),
      'THREE': cfvo.select('cfvo_30-60cm_mean'),
      'FOUR': cfvo.select('cfvo_60-100cm_mean'),
      'FIVE': cfvo.select('cfvo_100-200cm_mean')
}).clip(India_boundary).divide(10).rename('cfvo_f')

clay_f = clay.expression(
    '(0.005 * (0.5*((   (5*(ZERO+ONE)) + (10*(ONE+TWO)) + (15*(TWO + THREE))  +  (30*(THREE + FOUR))  +  (40*(FOUR + FIVE))  +   (100*(FIVE))    ) )  )) ', {
      'ZERO': clay.select('clay_0-5cm_mean'),
      'ONE': clay.select('clay_5-15cm_mean'),
      'TWO': clay.select('clay_15-30cm_mean'),
      'THREE': clay.select('clay_30-60cm_mean'),
      'FOUR': clay.select('clay_60-100cm_mean'),
      'FIVE': clay.select('clay_100-200cm_mean')
}).clip(India_boundary).divide(10).rename('clay_f')

nitrogen_f = nitrogen.expression(
    '(0.005 * (0.5*((   (5*(ZERO+ONE)) + (10*(ONE+TWO)) + (15*(TWO + THREE))  +  (30*(THREE + FOUR))  +  (40*(FOUR + FIVE))  +   (100*(FIVE))    ) )  )) ', {
      'ZERO': nitrogen.select('nitrogen_0-5cm_mean'),
      'ONE': nitrogen.select('nitrogen_5-15cm_mean'),
      'TWO': nitrogen.select('nitrogen_15-30cm_mean'),
      'THREE': nitrogen.select('nitrogen_30-60cm_mean'),
      'FOUR': nitrogen.select('nitrogen_60-100cm_mean'),
      'FIVE': nitrogen.select('nitrogen_100-200cm_mean')
}).clip(India_boundary).divide(100).rename('nitrogen_f')

phh2o_f = phh2o.expression(
    '(0.005 * (0.5*((   (5*(ZERO+ONE)) + (10*(ONE+TWO)) + (15*(TWO + THREE))  +  (30*(THREE + FOUR))  +  (40*(FOUR + FIVE))  +   (100*(FIVE))    ) )  )) ', {
      'ZERO': phh2o.select('phh2o_0-5cm_mean'),
      'ONE': phh2o.select('phh2o_5-15cm_mean'),
      'TWO': phh2o.select('phh2o_15-30cm_mean'),
      'THREE': phh2o.select('phh2o_30-60cm_mean'),
      'FOUR': phh2o.select('phh2o_60-100cm_mean'),
      'FIVE': phh2o.select('phh2o_100-200cm_mean')
}).clip(India_boundary).divide(10).rename('phh2o_f')

silt_f = silt.expression(
    '(0.005 * (0.5*((   (5*(ZERO+ONE)) + (10*(ONE+TWO)) + (15*(TWO + THREE))  +  (30*(THREE + FOUR))  +  (40*(FOUR + FIVE))  +   (100*(FIVE))    ) )  )) ', {
      'ZERO': silt.select('silt_0-5cm_mean'),
      'ONE': silt.select('silt_5-15cm_mean'),
      'TWO': silt.select('silt_15-30cm_mean'),
      'THREE': silt.select('silt_30-60cm_mean'),
      'FOUR': silt.select('silt_60-100cm_mean'),
      'FIVE': silt.select('silt_100-200cm_mean')
}).clip(India_boundary).divide(10).rename('silt_f')

sand_f = sand.expression(
    '(0.005 * (0.5*((   (5*(ZERO+ONE)) + (10*(ONE+TWO)) + (15*(TWO + THREE))  +  (30*(THREE + FOUR))  +  (40*(FOUR + FIVE))  +   (100*(FIVE))    ) )  )) ', {
      'ZERO': sand.select('sand_0-5cm_mean'),
      'ONE': sand.select('sand_5-15cm_mean'),
      'TWO': sand.select('sand_15-30cm_mean'),
      'THREE': sand.select('sand_30-60cm_mean'),
      'FOUR': sand.select('sand_60-100cm_mean'),
      'FIVE': sand.select('sand_100-200cm_mean')
}).clip(India_boundary).divide(10).rename('sand')

soc_f = soc.expression(
    '(0.005 * (0.5*((   (5*(ZERO+ONE)) + (10*(ONE+TWO)) + (15*(TWO + THREE))  +  (30*(THREE + FOUR))  +  (40*(FOUR + FIVE))  +   (100*(FIVE))    ) )  )) ', {
      'ZERO': soc.select('soc_0-5cm_mean'),
      'ONE': soc.select('soc_5-15cm_mean'),
      'TWO': soc.select('soc_15-30cm_mean'),
      'THREE': soc.select('soc_30-60cm_mean'),
      'FOUR': soc.select('soc_60-100cm_mean'),
      'FIVE': soc.select('soc_100-200cm_mean')
}).clip(India_boundary).divide(10).rename('soc')


ocd_f = ocd.expression(
    '(0.005 * (0.5*((   (5*(ZERO+ONE)) + (10*(ONE+TWO)) + (15*(TWO + THREE))  +  (30*(THREE + FOUR))  +  (40*(FOUR + FIVE))  +   (100*(FIVE))    ) )  )) ', {
      'ZERO': ocd.select('ocd_0-5cm_mean'),
      'ONE': ocd.select('ocd_5-15cm_mean'),
      'TWO': ocd.select('ocd_15-30cm_mean'),
      'THREE': ocd.select('ocd_30-60cm_mean'),
      'FOUR': ocd.select('ocd_60-100cm_mean'),
      'FIVE': ocd.select('ocd_100-200cm_mean')
}).clip(India_boundary).divide(10).rename('ocd')

# Exporting data to drive
geemap.ee_export_image_to_drive(
   bdod_f , 
   description="bdod_f_india", 
   folder="gee",scale=10 , 
   region=India_boundary.geometry().dissolve() , 
   maxPixels=1e13
)


# Processing for climate data extraction
AET_A = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate('2023-01-01','2023-12-30').select('aet').mean().multiply(0.1).clip(India_boundary)
DEF_A = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate('2023-01-01','2023-12-30').select('def').mean().multiply(0.1).clip(India_boundary)
PDSI_A = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate('2023-01-01','2023-12-30').select('pdsi').mean().multiply(0.01).clip(India_boundary)
PR_A = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate('2023-01-01','2023-12-30').select('pr').mean().clip(India_boundary)
SOIL_A = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate('2023-01-01','2023-12-30').select('soil').mean().multiply(0.01).clip(India_boundary)
SRAD_A = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate('2023-01-01','2023-12-30').select('srad').mean().multiply(0.01).clip(India_boundary)
TMIN_A = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate('2023-01-01','2023-12-30').select('tmmn').mean().multiply(0.01).clip(India_boundary)
TMAX_A = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate('2023-01-01','2023-12-30').select('tmmx').mean().multiply(0.01).clip(India_boundary)

# Exporting data to drive
geemap.ee_export_image_to_drive(
   TMAX_A , 
   description="TMAX_A_india", 
   folder="gee",
   scale=10 , 
   region=India_boundary.geometry().dissolve() , 
   maxPixels=1e13
)

# Processing for Topography data extraction
elevation = ee.Image('USGS/SRTMGL1_003').select('elevation').rename('elevation').updateMask(masked_lc)
slope = ee.Terrain.slope(elevation).rename('slope').rename('slope').updateMask(masked_lc)
aspect = ee.Terrain.aspect(elevation).rename('aspect').rename('aspect').updateMask(masked_lc)

# Exporting data to drive
geemap.ee_export_image_to_drive(
   aspect , 
   description="aspect_india", 
   folder="gee",
   scale=10 , 
   region=India_boundary.geometry().dissolve() , 
   maxPixels=1e13
)






# Model 

import geopandas as gpd
import rasterio
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, GridSearchCV
from rasterio.windows import Window
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm  # Import tqdm for progress bars

# Function to extract raster values at point locations in chunks with progress bar
def extract_raster_values_chunked(raster_paths, points, chunk_size):
    all_values = []
    valid_indices = []

    # Adding a progress bar for the point iteration
    for idx, point in tqdm(enumerate(points.geometry), total=len(points), desc="Extracting raster values"):
        point_values = []
        valid_point = True
        for raster_path in raster_paths:
            with rasterio.open(raster_path) as src:
                try:
                    row, col = src.index(point.x, point.y)
                    # Define a small window around the point
                    row_start = max(0, row - chunk_size // 2)
                    col_start = max(0, col - chunk_size // 2)
                    window = Window(col_start, row_start, chunk_size, chunk_size)
                    data = src.read(1, window=window)

                    value = data[row - row_start, col - col_start]
                    # Check if the pixel is within valid range (not no-data)
                    if src.nodata is not None and value == src.nodata:
                        valid_point = False
                    point_values.append(value)
                except IndexError:
                    valid_point = False  # Exclude points outside raster extent

        if valid_point:
            all_values.append(point_values)
            valid_indices.append(idx)

    return np.array(all_values), valid_indices

# Function to predict and save a new raster in chunks with progress bar
def predict_new_raster_chunked(model, new_raster_paths, output_path, chunk_size, feature_columns):
    with rasterio.open(new_raster_paths[0]) as src:
        profile = src.profile
        profile.update(dtype=rasterio.float64, count=1)  # Use float64 for larger NoData value

        rows, cols = src.shape
        nodata_value = src.nodata  # Get the no-data value from the raster


      # Counter to track the number of DataFrames printed
    printed_chunks = 0
    max_chunks_to_print = 3  # Number of chunk DataFrames to print

    # Create a progress bar for chunk processing
    with rasterio.open(output_path, 'w', **profile) as dst:
        total_chunks = (rows // chunk_size + 1) * (cols // chunk_size + 1)
        with tqdm(total=total_chunks, desc="Predicting in chunks") as pbar:
            for row_start in range(0, rows, chunk_size):
                for col_start in range(0, cols, chunk_size):
                    row_end = min(row_start + chunk_size, rows)
                    col_end = min(col_start + chunk_size, cols)

                    # Prepare a window for the current chunk
                    window = Window(col_start, row_start, col_end - col_start, row_end - row_start)
                    chunk_features = []

                    for raster_path in new_raster_paths:
                        with rasterio.open(raster_path) as src:
                            data = src.read(1, window=window)

                            # Mask the no-data value
                            if nodata_value is not None:
                                data[data == nodata_value] = np.nan  # Set no-data pixels to NaN

                            chunk_features.append(data.flatten())

                    chunk_features = np.array(chunk_features).T
                    predictions = np.full(chunk_features.shape[0], np.nan, dtype=np.float64)  # Use float64 for predictions

                    # Convert chunk features to a DataFrame with feature names
                    chunk_features_df = pd.DataFrame(chunk_features, columns=feature_columns)

                  

                    # Identify valid pixels (pixels that are not NaN)
                    valid_pixels = np.isfinite(chunk_features_df).all(axis=1)  # Exclude NaN values (no-data pixels)


  # Print the chunk only if it has valid values
                    if np.any(valid_pixels) and printed_chunks < max_chunks_to_print:
                        print(f"Chunk DataFrame #{printed_chunks + 1} with valid values:")
                        print(chunk_features_df[valid_pixels].head())  # Print only valid rows
                        printed_chunks += 1

                    
                    if np.any(valid_pixels):
                        predictions[valid_pixels] = model.predict(chunk_features_df[valid_pixels])

                    # Reshape predictions to match the chunk size and write to output
                    predicted_chunk = predictions.reshape((row_end - row_start, col_end - col_start))
                    dst.write(predicted_chunk, 1, window=window)

                    pbar.update(1)  # Update progress bar



# Load training data 
shapefile_path = "path\\tree.shp"
raster_paths = ["raster1.tif", "raster1.tiff"]
points = gpd.read_file(shapefile_path)
class_values = points['classvalue']

# Select 100 points per class
points_per_class = points

# Extract raster values in chunks with progress bar
chunk_size = 512  # Define the chunk size (adjust as needed)
features, valid_indices = extract_raster_values_chunked(raster_paths, points_per_class, chunk_size)
filtered_classes = points_per_class['classvalue'].iloc[valid_indices]

# Include additional columns from the shapefile (Lat, Long, classname, etc.)
additional_columns = points_per_class[['Lat', 'Long', 'classname']].iloc[valid_indices].reset_index(drop=True)

# Min-Max normalization of the features
scaler = MinMaxScaler()
features = scaler.fit_transform(features)  # Min-Max scale the raster features

# Create DataFrame for model training
data = pd.DataFrame(features, columns=[f"Feature_{i}" for i in range(len(raster_paths))])
data['Class'] = filtered_classes.values
data['Lat'] = additional_columns['Lat']
data['Long'] = additional_columns['Long']
data['Classname'] = additional_columns['classname']

# Train-test split
X = data.drop(columns=["Class",'Lat','Long','Classname'])
y = data["Class"]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Hyperparameter tuning for Random Forest with progress bar
param_grid = {
    'n_estimators': [100, 200, 500],
    'max_depth': [10, 20, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'bootstrap': [True, False]
}
rf = RandomForestClassifier(random_state=42)
grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=3, n_jobs=-1, verbose=2)

# Adding a progress bar for the training process
with tqdm(total=1, desc="Hyperparameter tuning") as pbar:
    grid_search.fit(X_train, y_train)
    pbar.update(1)

# Best model
best_model = grid_search.best_estimator_
print(f"Best Parameters: {grid_search.best_params_}")
print(f"Model Accuracy: {best_model.score(X_test, y_test):.2f}")


from sklearn.metrics import classification_report, confusion_matrix, accuracy_score

# Best model
best_model = grid_search.best_estimator_
print(f"Best Parameters: {grid_search.best_params_}")

# Predict the test set
y_pred = best_model.predict(X_test)

# Overall Accuracy
overall_accuracy = accuracy_score(y_test, y_pred)
print(f"Overall Model Accuracy: {overall_accuracy:.2f}")

# Confusion Matrix
conf_matrix = confusion_matrix(y_test, y_pred)
print("Confusion Matrix:")
print(conf_matrix)

# Detailed Classification Report
class_names = list(map(str, y.unique()))  # Ensure class names are strings
class_report = classification_report(y_test, y_pred, target_names=class_names)
print("Classification Report:")
print(class_report)

# Class-Wise Accuracy
class_accuracies = conf_matrix.diagonal() / conf_matrix.sum(axis=1)  # Diagonal/row sums
print("\nClass-Wise Accuracy:")
for class_name, accuracy in zip(class_names, class_accuracies):
    print(f"Class {class_name} Accuracy: {accuracy:.2f}")


# Predict new raster in chunks with progress bar
new_raster_paths = ["new_raster1.tif", "new_raster2.tif"]
output_path = "path\\prediction.tif"

feature_columns = [f"Feature_{i}" for i in range(len(raster_paths))] 
predict_new_raster_chunked(best_model, new_raster_paths, output_path, chunk_size, feature_columns)













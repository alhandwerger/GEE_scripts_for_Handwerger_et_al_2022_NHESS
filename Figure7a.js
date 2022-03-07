/*
MIT License

Copyright (c) 2021 Mong-Han Huang and Alexander L. Handwerger

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
// #############################################################################
// ### Event info ###
// #############################################################################
/* 
/* 
* SAR-based landslide (and other ground surface) detection 
* Study Site: Huong Phung, Vietnam
* Event: high rainfall 18 October 2020
* Code written by Mong-Han Huang, Department of Geology, mhhuang@umd.edu
* and Alexander L. Handwerger, UCLA JIFRESSE/NASA JPL, alhandwerger@g.ucla.edu and alexander.handwerger@jpl.nasa.gov with help from many users on stackexchange!
*/   
 
// #############################################################################
// ### Set AOI ###
// #############################################################################
var AOI = AOI; // define Area of Interest
Map.setOptions("TERRAIN");
Map.centerObject(AOI, 12); //zooms to center of AOI after clicking "run". The number determines the zoom level.
 
// #############################################################################
// ### Define pre-event and post-event time periods ###
// #############################################################################

//**************** Define pre-event and post-event time periods ********************************//
//define pre-event stack time period
var PreEventTime_1 = '2016-10-01T23:59'; // format: yyyy-mm-dd-HH:MM
var PreEventTime_2 = '2020-10-17T23:59'; // format: yyyy-mm-dd-HH:MM

//define post-event stack time period
var PostEventTime_1 = '2020-10-19T23:59'; // format: yyyy-mm-dd-HH:MM
var PostEventTime_2 = '2020-11-01T23:59'; // format: yyyy-mm-dd-HH:MM
var PostEventTime_S2longer = '2021-10-25T23:59'; // format: yyyy-mm-dd-HH:MM

// #############################################################################
// ### slope and curvature thresholds ###
// #############################################################################

// DEM thresholds
var slope_threshold = 12;    // Exclude areas with hillslope angle < slope_threshold, unit: degree
var curv_threshold = -0.005; // Exclude areas with hillslope curvature > curv_threshold, unit: m/m^2


// #############################################################################
// ### NASA DEM ###
// #############################################################################


// load NASA Digital Elevation Model (DEM)
var NASADEM_dataset = ee.Image('NASA/NASADEM_HGT/001')

var elevation = NASADEM_dataset.select('elevation');
elevation=elevation.clip(AOI); //clip to AOI

var waterZones = NASADEM_dataset.select('swb'); //load in water body data
var waterMask = waterZones.eq(0) // Create a binary water mask.
var exaggeration = 1; //changes the appearance of the hillshade
var hillshade = ee.Terrain.hillshade(elevation.multiply(exaggeration));
hillshade=hillshade.clip(AOI);
var slope = ee.Terrain.slope(elevation);          // slope angle in degrees
slope=slope.clip(AOI);
var mask_slope = slope.gte(slope_threshold);      // slope mask with values 0 or 1, removes values less than or equal to threshold
var slope_masked = slope.updateMask(mask_slope);  // slope angle with values < threshold excluded

// Calculate hillslope curvature
// Define a Gaussian kernel for smoothing. This step helps reduce noise in the curvature maps
var smooth_curv = ee.Kernel.gaussian({
  radius: 60,
  sigma: 30,
  units: 'meters',
  normalize: true,
});

// Smoothing the DEM with the gaussian kernel.
var elevation_smooth= elevation.convolve(smooth_curv).resample("bilinear");
var xyDemGrad = elevation_smooth.gradient().resample("bilinear");
var xGradient = xyDemGrad.select('x').gradient().resample("bilinear");
var yGradient = xyDemGrad.select('y').gradient().resample("bilinear");
var curvature = xGradient.select('x').add(yGradient.select('y'));
var mask_curvature = curvature.gte(curv_threshold);
var curvature_masked = curvature.updateMask(mask_curvature);
//var DEMMask= mask_slope.multiply(mask_curvature).multiply(waterMask); //slope, curvature, and water mask

// #############################################################################
// ### Sentinel-2 Optical Images ###
// #############################################################################

// cloud filter and mask for Sentinel-2 optical data
function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}
// Load Sentinel-2 reflectance data.
var S2_PreEvent = ee.ImageCollection('COPERNICUS/S2')
                  .filterDate(PreEventTime_1,PreEventTime_2)
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) //only include images with less than 10% clouds
                  .filterBounds(AOI)
                  .map(maskS2clouds);
                  
var S2_PreEvent_median=S2_PreEvent.median();
S2_PreEvent_median=S2_PreEvent_median.clip(AOI);

var S2_PostEvent = ee.ImageCollection('COPERNICUS/S2')
                  .filterDate(PostEventTime_1, PostEventTime_S2longer)
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) //only include images with less than 10% 
                  .filterBounds(AOI)
                  .map(maskS2clouds); 
                  
var S2_PostEvent_median = S2_PostEvent.median();
S2_PostEvent_median=S2_PostEvent_median.clip(AOI);

// #############################################################################
// ### Sentinel-1 SAR Images ###
// #############################################################################

// LOAD Sentinel-1 (S1) Backscatter Coefficient data VH polarization
var imgVH = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .select('VH') // VH polarization only
        .filterBounds(AOI)
        .map(function(image) {
          var edge = image.lt(-30.0); //remove low edge values as suggested by GEE
          var maskedImage = image.mask().and(edge.not());
          return image.updateMask(maskedImage);
        });

var desc = imgVH.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')); //descending acquisition geometry data
var asc = imgVH.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));  //ascending acquisition geometry data

var PreEventPeriod = ee.Filter.date(PreEventTime_1,PreEventTime_2);
var PostEventPeriod = ee.Filter.date(PostEventTime_1,PostEventTime_2);

// calculate the median of the Pre-event S1 SAR Backscatter Intensity
var PreEventPeriod_asc = ee.Image.cat(asc.filter(PreEventPeriod).median());
var PreEventPeriod_desc = ee.Image.cat(desc.filter(PreEventPeriod).median());

// calculate the median of the Post-event S1 SAR Backscatter Intensity
var PostEventPeriod_asc = ee.Image.cat(asc.filter(PostEventPeriod).median());
var PostEventPeriod_desc = ee.Image.cat(desc.filter(PostEventPeriod).median());

// print out image information (number of images, image file name) 
//Note this counts two images acquired on the same day as two separate images, to avoid double counting create very small AOI called "AOI_4_totalimagecount"
var num_asc_pre = asc.filter(PreEventPeriod).filterBounds(AOI_4_totalimagecount);
var num_desc_pre = desc.filter(PreEventPeriod).filterBounds(AOI_4_totalimagecount);
var num_asc_post = asc.filter(PostEventPeriod).filterBounds(AOI_4_totalimagecount);
var num_desc_post = desc.filter(PostEventPeriod).filterBounds(AOI_4_totalimagecount);
var count_asc_pre = num_asc_pre.sort('system:time_start').toList(5000,0).length();   // 5000 controls size of the list
var count_desc_pre = num_desc_pre.sort('system:time_start').toList(5000,0).length(); // 
var count_asc_post = num_asc_post.sort('system:time_start').toList(5000,0).length(); // 
var count_desc_post = num_desc_post.sort('system:time_start').toList(5000,0).length(); // 

print("Number of pre-event ascending images: ", count_asc_pre  );
print("Number of pre-event descending images: ", count_desc_pre  );
print("Number of post-event ascending images: ", count_asc_post  );
print("Number of post-event descending images: ", count_desc_post  );
print(num_asc_post,'Post-event ascending');   //print out to check date of image collection
print(num_desc_post,'Post-event descending'); //print out to check date of image collection


// #############################################################################
// ### I_ratio ###
// #############################################################################

// calculate the log ratio (using subtraction since data are in log scale) for Pre- and Post-event S1 SAR Backscatter
var I_ratio_desc = PreEventPeriod_desc.subtract(PostEventPeriod_desc);
I_ratio_desc=I_ratio_desc.clip(AOI);

var I_ratio_asc = PreEventPeriod_asc.subtract(PostEventPeriod_asc);
I_ratio_asc=I_ratio_asc.clip(AOI);

var I_ratio_avg_desc_asc = (I_ratio_asc.add(I_ratio_desc)).divide(2); // calculate the mean backscatter change for ascending and descending scenes combined


// #############################################################################
// ### Define the I_ratio to be used for landslide detection. This can be "I_ratio_desc", "I_ratio_asc", or "I_ratio_avg_desc_asc". 
// ### We recommend using "I_ratio_avg_desc_asc" if both ascending and descending data exist. Otherwise use either "I_ratio_desc" or "I_ratio_asc" ###
// #############################################################################

//slope and water mask
//var I_ratio_avg_masked = I_ratio_avg_desc_asc.updateMask(mask_slope).updateMask(waterMask); 

//slope, curvature, and water mask
var I_ratio_avg_masked = I_ratio_avg_desc_asc.updateMask(mask_slope).updateMask(mask_curvature).updateMask(waterMask); 

// #############################################################################
// ### Calculate I_ratio statistics ###
// #############################################################################

// Calculate percentiles of I_ratio
var I_ratio_Percentiles = I_ratio_avg_masked.reduceRegion({
  reducer:ee.Reducer.percentile([90, 95, 99]),
  geometry: AOI,
  scale: 10,
  bestEffort: true
});
print('I_ratio Percentiles',I_ratio_Percentiles)

var I_ratio_90thPtile = I_ratio_avg_masked.gte(ee.Number(I_ratio_Percentiles.get("VH_p90")))
var I_ratio_95thPtile = I_ratio_avg_masked.gte(ee.Number(I_ratio_Percentiles.get("VH_p95")))
var I_ratio_99thPtile = I_ratio_avg_masked.gte(ee.Number(I_ratio_Percentiles.get("VH_p99")))


// #############################################################################
// ### Heatmap and threshold-based landslide detection ###
// #############################################################################

var backscatter_chng_threshold = ee.Number(I_ratio_Percentiles.get("VH_p99")); // backscatter change threshold. This determines the values that correspond to potential landslides. Used for heatmap, polygons, points, and raster of possible landslide locations
print('I_ratio threshold',I_ratio_Percentiles.get("VH_p99"))

//  identify possible landslide zones
var pos_slide_zones = I_ratio_avg_masked.gt(backscatter_chng_threshold);
pos_slide_zones = pos_slide_zones.updateMask(pos_slide_zones.neq(0));

// Convert the pos_slide_zones of the thresholded to vectors.
var pos_slide_polygons = pos_slide_zones.addBands(I_ratio_avg_masked).reduceToVectors({
  geometry: AOI,
  crs: I_ratio_avg_masked.projection(),
  scale: 10,
  bestEffort: true,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'zone',
  reducer: ee.Reducer.mean()
});

// Convert from polygons to points
var S1_points =  imgVH;
var first = ee.Image(S1_points.first()).clip(pos_slide_polygons)
// get image projection
var proj = first.select([0]).projection()
// get coordinates image
var latlon = ee.Image.pixelLonLat().reproject(proj)
// put each lon lat in a list
var coords = latlon.select(['longitude', 'latitude'])
                 .reduceRegion({
  reducer: ee.Reducer.toList(),
  geometry: pos_slide_polygons,
  scale: 10,
  bestEffort: true,
})

// get lat & lon
var lat = ee.List(coords.get('latitude'))
var lon = ee.List(coords.get('longitude'))

// zip them. 
var point_list = lon.zip(lat)

// Create points
var mp = ee.Geometry.MultiPoint(point_list)

//******************** Identify points for heatmap to be exported as a KML ***********************//

var mp2 = ee.FeatureCollection(point_list.map(function(p){
  var point = ee.Feature(ee.Geometry.Point(p), {})
  return point
}))

var points_heatmap = mp2.map(function(feature){
  return feature.set('dummy',1);
});


// #############################################################################
// ### Color Palettes for Maps ###
// #############################################################################
// define color palette ('ffffff' is white, 'ff0000' is red, '0000ff' is blue, '000000' is black)
var ColorScale = {min: -10, max: 10, palette: ['0013ff','8178ff','ffffff','ff7e7e','ff0000']}; // for Backscatter Intensity change (I_ratio)
var rgbVis = {min: 0.0, max: 0.18, bands: ['B4', 'B3', 'B2']}; // for sentinel-2
var rgbVis_wClouds = {min: 0.0, max: 2000, bands: ['B4', 'B3', 'B2']}; // for sentinel-2 with clouds
var ColorCurv = {min: -0.02, max:0.02, palette: ['ff0000','ffffff','0000ff']}; // for curvature 
var MaskColor = {min: 0, max:1, palette: ['ffffff','0050d6']}; // for curvature 
var percentileColor = {min: 0, max:1, palette: ['ffffff','c000d6']}; // for showing percentile  

// #############################################################################
// ### Add layers to maps ###
// #############################################################################

//DEM layers
Map.addLayer(hillshade, null, 'NASADEM Hillshade');
//topo lines
var lines = ee.List.sequence(0, 4000, 100)
//var lines = ee.List.sequence(0, 4000, 200)
var contourlines = lines.map(function(line) {
  var mycontour = elevation
    .convolve(ee.Kernel.gaussian(5, 3))
    .subtract(ee.Image.constant(line)).zeroCrossing() 
    .multiply(ee.Image.constant(line)).toFloat();
  return mycontour.mask(mycontour);
})
contourlines = ee.ImageCollection(contourlines).mosaic()
Map.addLayer(contourlines, {min: 0, max: 1000, palette:['000000', '000000']}, 'topo lines', false, 0.5)

//Sentinel-2 optical images
Map.addLayer(S2_PreEvent_median, rgbVis , 'S2 pre-event',false);
Map.addLayer(S2_PostEvent_median, rgbVis, 'S2 post-event',false);

//Possible landslide areas defined based on backscatter change threshold
Map.addLayer(pos_slide_zones.clip(AOI), {min: 0, max: 1, palette: ['ffffff','FF0000']}, 'possible landslide zones',false);
Map.addLayer(I_ratio_avg_masked, ColorScale, 'I_ratio masked', false);
//Map.addLayer(I_ratio_90thPtile, percentileColor, 'I_ratio >= 90th percentile', false);
Map.addLayer(I_ratio_95thPtile, percentileColor, 'I_ratio >= 95th percentile', false);
Map.addLayer(I_ratio_99thPtile, percentileColor, 'I_ratio >= 99th percentile', false, 0.75);

// #############################################################################
// ### Generate heatmap visualization in GEE ###
// #############################################################################

//NOTE This often causes GEE to timeout. Until we can optimize the code, we leave this commented out and perform the heatmap calculation in QGIS ***********************//

var radius_heatmap2 = 20;         // heatmap kernal radius

function heatmap(mp2,radius_heatmap){
  var ptImg = points_heatmap.reduceToImage(['dummy'],ee.Reducer.first()).unmask(0);
  var kernel = ee.Kernel.circle(radius_heatmap);
  var result = ptImg.convolve(kernel);
  return result.clip(AOI); //clip to AOI
  //return result.updateMask(result.neq(0)); //uncomment this if want to remove green area
}

var heatmapImg = heatmap(points_heatmap,radius_heatmap2); 

var gradientcolormap = ['lightgreen','yellow','red']; //heatmap gradient
var heatmapColor = {
  min: 0,
  max: .1,
  palette: gradientcolormap,
};
//heatmap
Map.addLayer(heatmapImg,heatmapColor,'Heatmap',true, 0.6);

// #############################################################################
// ### save data for export ###
// #############################################################################

Export.image.toDrive({
  image: S2_PostEvent_median,
  description: 'HuongPhung_S2_PostEvent',
  scale: 10,
  fileFormat: 'GeoTIFF',
  region: AOI
});

Export.image.toDrive({
  image: S2_PreEvent_median,
  description: 'HuongPhung_S2_PreEvent',
  scale: 10,
  fileFormat: 'GeoTIFF',
  region: AOI
});

Export.image.toDrive({
  image: I_ratio_99thPtile,
  description: 'HuongPhung_I_ratio>=99thPtile',
  scale: 10,
  fileFormat: 'GeoTIFF',
  region: AOI
});

Export.image.toDrive({
  image: I_ratio_avg_masked,
  description: 'HuongPhung_I_ratio_masked',
  scale: 10,
  fileFormat: 'GeoTIFF',
  region: AOI
});

Export.table.toDrive({
  collection: points_heatmap,
  description:'HuongPhung_Points4Heatmap',
  fileFormat: 'KML'
});




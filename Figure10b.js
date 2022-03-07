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
* Study Site: Hokkaido, Japan 
* Event: no major landslide event  
* Code written by Mong-Han Huang, Department of Geology, mhhuang@umd.edu
* and Alexander L. Handwerger, UCLA JIFRESSE/NASA JPL, alhandwerger@g.ucla.edu and alexander.handwerger@jpl.nasa.gov with help from many users on stackexchange!
*/   
  
// #############################################################################
// ### Set AOI ###
// #############################################################################
var AOI = Sub_AOI; // define Area of Interest
Map.setOptions("TERRAIN");
Map.centerObject(AOI, 15); //zooms to center of AOI after clicking "run". The number determines the zoom level.
 
// #############################################################################
// ### Define pre-event and post-event time periods ### 
// #############################################################################
//define pre-event stack time period
var PreEventTime_1 = '2017-11-01T23:59'; // format: yyyy-mm-dd-HH:MM
var PreEventTime_2 = '2018-06-05T23:59'; // format: yyyy-mm-dd-HH:MM

//define post-event stack time period for event inventory detection
//for this approach use all of the available post-event imagery
var PostEventTime_1 = '2018-06-06T23:59'; // format: yyyy-mm-dd-HH:MM
var PostEventTime_2 = '2018-06-21T23:59'; // format: yyyy-mm-dd-HH:MM
var PostEventTime_S2longer = '2018-07-25T23:59'; // format: yyyy-mm-dd-HH:MM


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
// ### Color Palettes for Maps ###
// #############################################################################
// define color palette ('ffffff' is white, 'ff0000' is red, '0000ff' is blue, '000000' is black)
var rgbVis = {min: 0.0, max: 0.18, bands: ['B4', 'B3', 'B2']}; // for sentinel-2
var rgbVis_1 = {min: 0.0, max: 0.3, bands: ['B4', 'B3', 'B2']}; // for sentinel-2 with clouds

// #############################################################################
// ### Add layers to maps ###
// #############################################################################

//Sentinel-2 optical images
Map.addLayer(S2_PreEvent_median, rgbVis_1 , 'S2 pre-event',true);
Map.addLayer(S2_PostEvent_median, rgbVis, 'S2 post-event',false);


// #############################################################################
// ### save data for export ###
// #############################################################################

Export.image.toDrive({
  image: S2_PostEvent_median,
  description: 'Hokkaido_S2_PostEvent',
  scale: 10,
  fileFormat: 'GeoTIFF',
  region: AOI
});

Export.image.toDrive({
  image: S2_PreEvent_median,
  description: 'Hokkaido_S2_PreEvent',
  scale: 10,
  fileFormat: 'GeoTIFF',
  region: AOI
});

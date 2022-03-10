// The following code estimates projected habitat-suitable range (HSR) loss to urban and non-urban land-use change between 2015 and 2050  as per the analysis undertaken in Simkin et al. (2022):
//
// Simkin, R.D., Seto, K.S., McDonald, R.I., Jetz, W., 2022, Biodiversity impacts and conservation implications of urban land expansion projected to 2050,
// Proc. Natl. Acad. Sci. U.S.A. (X), XX-XX, {DOI}
//
//
// The code provided below estimates HSR loss for a single species but can be adapted to larger groups as required.
// For full details of the method refer to the publication.

// Questions can be directed to rohan.simkin@yale.edu



/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//IMPORT DATA
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

// Countries
// Natural Earth 50m countries dataset www.naturalearthdata.com
var countries = ee.FeatureCollection('projects/map-of-life/habitat_projection/urban/countries');

// Cities
// Natural Earth 10m populated places dataset www.naturalearthdata.com
var cities = ee.FeatureCollection('projects/map-of-life/habitat_projection/urban/cities');


// GLOBIO LULC images
// Available at https://www.globio.info/globio-data-downloads 
var globio_2015 = ee.Image('projects/map-of-life/habitat_projection/urban/globio_2015');
var globio_ssp1_2050 = ee.Image('projects/map-of-life/habitat_projection/urban/globio_ssp1_2050');
// Combined into a single image with bands '2015' and '2050'
var globio_ssp1 = ee.ImageCollection.fromImages([globio_2015, globio_ssp1_2050])
  .toBands()
  .select(['0_b1', '1_b1'], ['2015', '2050']);

var globio_ssp3_2050 = ee.Image('projects/map-of-life/habitat_projection/urban/globio_ssp3_2050');
var globio_ssp3 = ee.ImageCollection.fromImages([globio_2015, globio_ssp3_2050])
  .toBands()
  .select(['0_b1', '1_b1'], ['2015', '2050']);

var globio_ssp5_2050 = ee.Image('projects/map-of-life/habitat_projection/urban/globio_ssp5_2050');
var globio_ssp5 = ee.ImageCollection.fromImages([globio_2015, globio_ssp5_2050])
  .toBands()
  .select(['0_b1', '1_b1'], ['2015', '2050']);
  
// Set the scale for any reductions undertaken
// This is the native scale of the GLOBIO LULC dataset
var scale = ee.Number(globio_2015.projection().nominalScale())


// URBAN LAND
// Import urban forecasts
var urban_2015 = ee.Image('projects/map-of-life/habitat_projection/urban/urban-2015').unmask(0);
var ssp1 = ee.Image('projects/map-of-life/habitat_projection/urban/urban-ssp1').unmask(0);
var ssp3 = ee.Image('projects/map-of-life/habitat_projection/urban/urban-ssp3').unmask(0);
var ssp5 = ee.Image('projects/map-of-life/habitat_projection/urban/urban-ssp5').unmask(0);
Map.addLayer(ssp5.selfMask(), {palette: ['FCECF8', 'FC01BD']}, 'urbanLand_ssp5_2050', false);

// The DEM
var DEM = ee.Image("projects/map-of-life/habitat_projection/urban/DEM");

//The GLOBIO LULC forecasts are based upon the ESA's CCI Landcover dataset for 2015 - https://www.esa-landcover-cci.org/ 
//read in the CCI Landcover data
var cci_2015 = ee.Image('projects/map-of-life/habitat_projection/urban/cci_2015');

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//SPECIES RANGE AND HABITAT PREFERENCES
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

//Get the expert range map for the species
var sampleSpp = ee.ImageCollection(ee.Image("projects/map-of-life/habitat_projection/urban/demo/Oreophasis_derbianus"));

//Define habitat preferences for species of interest
var sampleSpp_prefs = ee.Dictionary({
elev_max: 3400, //Elevation range
elev_min: 1100, 
//Habitat preferences
habitats_lst: [1, //Evergreen Needleleaf Forests 
  2, //Evergreen Broadleaf Forests 
  3, //Deciduous Needleleaf Forests
  4, //Deciduous Broadleaf Forests 
  5] //Mixed Forests   
});

//Add the habitat prefs to the range map image
sampleSpp = sampleSpp.map(function(i) {return i.set(sampleSpp_prefs)});


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//RECLASSIFY PASTURE AND SECONDARY VEGETATION 
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

// This section creates a set of categorical LULC images that are used to identify HSR for species.

// For each species in our analysis we identified habitat preferences based on LULC classes defined by the International Geosphere-Biosphere Program (IGBP) land cover classification system consisting of 17 classes including natural land-covers and human land-uses
// Our assessment of habitat preferences did not include the LULC classes pasture, or secondary vegetation which are included in the GLOBIO LULC dataset.
// To account for this we match species to these LULC classes according to a set of assumptions outlined in Simkin et al. (2022).
// To identify species that may occur within low-intensity pasture and secondary vegetation pixels in the the GLOBIO LULC forecasts we first reclassify these pixels to equivilent LULC classes to which we can match species habitat preferences.
// High-intensity pasture pixels are retained at this stage and species are matched to (or excluded from) this LULC class at a later stage. 
 

// 2015 reclassified image
// For full details refer to Simkin et al. (2022)
// - Assigns all low-intensity pasture pixels back to the underlying CCI Landcover LULC category
// - All other pixels are unchanged

//A function that performs the reclassification
var pastureReclassifyFunction = function(img) {
  var pastureMask = img.neq(4); //mask: lowintensity pasture ==0, everything else ==1
  return img.updateMask(pastureMask).unmask(cci_2015); //replace pasture values with cci 2015 value
};
//The resulting reclassified 2015 GLOBIO LULC image
var globio_2015_pastureReclassified = pastureReclassifyFunction(globio_2015);

// 2050 reclassified images
// For full details refer to Simkin et al. (2022)
//  - Convert secondary vegetation pixels in 2050 back to the 2015 value (i.e. urban, crop, forestry or pasture)
//  - low intensity pasture in 2050 is reclassified based on its 2015 value
//      - Where it was natural in 2015, reclassify back to the natural land cover value
//      - where it was low intensity pasture in 2015, reclassify back to the underlying cci landcover
//      - where it was a non-natural land cover (crop, urban, high-intensity pasture) in 2015, retain as low intensity pasture and match species with habitat preference for IGBP LULC category 'Croplands' at later stage.

var noRegain_reclassifyFunction = function(img2015, img2050) {
  //Replace 2050 secondary vegetation pixels with the 2015 LULC class
  var secondaryMask = img2050.neq(6);
  var secondaryReclassified = img2050.mask(secondaryMask).unmask(img2015);
  //pasture pixel reclassification - takes the image with seconday pixels already reclassified as its input
  var reclassified = secondaryReclassified.where(
    secondaryReclassified.eq(4).and(img2015.gte(50).and(img2015.lt(230))), //where a pixel is both low intensity pasture in 2050 AND natural land cover in 2015
    img2015)//return the natural 2015 value
    .where(secondaryReclassified.eq(4).and(img2015.eq(4)), //where a pixel is low-intensity pasture in 2050 AND low-intensity pasture in 2015
    cci_2015);//return the cci natural land cover pixel value
  return reclassified; //The image with both pasture and secondary reclassified  
};

//Reclassify 2050 LULC images
var globio_ssp1_2050_noRegain = noRegain_reclassifyFunction(globio_2015, globio_ssp1_2050);
var globio_ssp3_2050_noRegain = noRegain_reclassifyFunction(globio_2015, globio_ssp3_2050);
var globio_ssp5_2050_noRegain = noRegain_reclassifyFunction(globio_2015, globio_ssp5_2050);


//Join the reclassified images for 2015 and 2050 into an image collection
// ssp1
var globio_ssp1_noRegain = ee.ImageCollection.fromImages([globio_2015_pastureReclassified, globio_ssp1_2050_noRegain])
  .toBands()
  .select(['0_b1', '1_b1'], ['2015', '2050']);

// ssp3
var globio_ssp3_noRegain = ee.ImageCollection.fromImages([globio_2015_pastureReclassified, globio_ssp3_2050_noRegain])
  .toBands()
  .select(['0_b1', '1_b1'], ['2015', '2050']);

// ssp5
var globio_ssp5_noRegain = ee.ImageCollection.fromImages([globio_2015_pastureReclassified, globio_ssp5_2050_noRegain])
  .toBands()
  .select(['0_b1', '1_b1'], ['2015', '2050']);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//URBAN CLUSTERS
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Identify clusters of urban land based on the forecasted extent of urban land in 2050
// Urban clusters are defined as:
// - contiguous groups of urban land pixels with probability of becoming urban > 0.25
// - must contain a city point from the Natural Earth Populated Places dataset
// - must be > 400 square km

//A function used to name the urban clusters based on the largest city within the cluster boundary
var cityNameFunction = function(feat) {
  var boundsFilter = ee.Filter.bounds(feat.geometry());//A filter that selects all cities within the feature
  var nCities = ee.Number(cities.filter(boundsFilter).size()); //no. cities that intersect the feature 
  var cityList = ee.List(ee.Algorithms.If({
    condition: nCities.lt(ee.Number(1)), //If no cities intersect
    trueCase: ee.FeatureCollection(ee.Feature(null).set('NAME', 'noCity').set('ADM0NAME', 'noCountry').set('ne_id', 'noId')).toList(1), //Return 'noCity'
    falseCase: ee.FeatureCollection(cities.filter(boundsFilter).sort({property: 'POP_MAX',ascending: false})).toList(10)
  }));//Get the 10 largest cities that intersect the cluster
  var cityNames = cityList.map(function(city){return ee.String(ee.Feature(city).get('NAME'))}).join(','); //Get city names as strings
  var cityIds = cityList.map(function(city){return ee.Number(ee.Feature(city).get('ne_id'))}).join(','); //Get city ids as numbers
  var countryOfLargestCity = ee.Feature(cityList.get(0)).getString('ADM0NAME'); //The country of the largest city within the cluster
  var ar = ee.Feature(feat).area(100); //Area of the cluster
  //Return a feature with properties set
  return ee.Feature(feat.geometry()).set('cities', cityNames)
    .set('countryOfLargestCity', countryOfLargestCity)
    .set('area', ar)
    .set('cityIds', cityIds);
};

//Identify and name urban clusters 
var ssp1Clusters = ssp1.gt(0.25).selfMask().reduceToVectors({
  geometry: ee.Geometry.Rectangle([-180,-85,180,85],'EPSG:4326',false),
  maxPixels: 1e12
}).map(cityNameFunction)
  .filter(ee.Filter.neq('cities', 'noCity'))
  .filter(ee.Filter.gt('area', 400000000));


var ssp3Clusters = ssp3.gt(0.25).selfMask().reduceToVectors({
  geometry: ee.Geometry.Rectangle([-180,-85,180,85],'EPSG:4326',false),
  maxPixels: 1e12
}).map(cityNameFunction)
  .filter(ee.Filter.neq('cities', 'noCity'))
  .filter(ee.Filter.gt('area', 400000000));


var ssp5Clusters = ssp5.gt(0.25).selfMask().reduceToVectors({
  geometry: ee.Geometry.Rectangle([-180,-85,180,85],'EPSG:4326',false),
  maxPixels: 1e12
}).map(cityNameFunction)
  .filter(ee.Filter.neq('cities', 'noCity'))
  .filter(ee.Filter.gt('area', 400000000));



/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//MERGE THE GLOBIO AND URBAN FORECASTS
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

// The GLOBIO LULC forecasts are used in our analysis to estimate non-urban HSR loss, while urban HSR losses are calculated seperately using dedicated urban land-use forecasts.
// As such, we create a version of the GLOBIO LULC forecast that holds urban land constant from 2015-2050

// In order to ensure that the area of urban land in the 2015 GLOBIO LULC image is consistent with the urban forecasts we merge the two datasets.
// The relassified GLOBIO LULC forecasts (desribed above) are used as inputs.
// The result is an image with bands for 2015 and 2050:
// - 2015 band contains the merged extent of urban land in the GLOBIO and urban forecasts for 2015.  All non urban land is unchanged from the 2015 reclassified GLOBIO LULC image.
// - Urban land extent in the 2050 band is unchanged from 2015.  Any land that changes from non-urban in 2015 to urban in 2050 in the GLOBIO image is reclassified back to its 2015 value. All non urban land is unchanged from the 2050 reclassified GLOBIO LULC image.

//Function that reclassifies the images
var globioUrbanMergeFunction = function(scenario) {
  //Create a mask that identifies combined urban land from urban and GLOBIO forecasts in 2015
  //Pixel values are 0 = urban, 1=non-urban  image is unmasked
  var urbanMask_2015 = (scenario.select('2015').eq(1).or(urban_2015.eq(1))).eq(0);
  
  //Create a version of the GLOBIO 2015 LULC image that masks out the merged urban extent for 2015
  //Image contains two bands (2015 and 2050), both masked to the extent of urban in 2015
  //i.e. the 2015 image contains no urban land, while the 2050 image contains only new urban land from the GLOBIO forecast
  var globio_2015UrbanMasked = scenario.updateMask(urbanMask_2015);
  
  //Create a version of the GLOBIO forecasts merged with the urban forecasts for 2015
  //The resulting image has bands for 2015 and 2050
  //Pixels in the 2050 image == 1 where they were urban in the 2015 GLOBIO OR urban images, plus any pixels that become urban in the globio image in 2050 
  var globio_urban2015Urban = globio_2015UrbanMasked.unmask(1);
  
  //Create a version of the urban globio merged layer that keeps urban area the same from 2015 - 2050
  //Any pixels in the 2050 band of the globio model that change from non-urban to urban are replaced with 
  //the non-urban pixel value that was present in 2015 - effectively keeping them constant
  var globio_urbanConstantUrban = globio_urban2015Urban.select('2050').updateMask(urbanMask_2015)//Get the 2050 image and mask out anything that was urban in 2015
    .eq(1).selfMask()//Identify just the pixels that were not urban in 2015 but become urban in 2050
    .multiply(globio_urban2015Urban.select('2015'))//Change those pixels back to thier 2015 value
    .unmask(globio_urban2015Urban.select('2050'))//Unmask all other pixels - they take on the 2050 globio value
    .addBands(globio_urban2015Urban.select('2015'));//Add the 2015 image as a band

  return globio_urbanConstantUrban;
};

// //Generate the merged datasets
// //Note that the 'base' scenario here refers to the images that do not have pasture and secondary vegetation reclassified.  These are used to determine the drivers of HSR change at a later stage.

//Base (unmodified) globio scenario
var ssp1_globio_constantUrban = globioUrbanMergeFunction(globio_ssp1);
//print('ssp1_globio_constantUrban', ssp1_globio_constantUrban)
var ssp3_globio_constantUrban = globioUrbanMergeFunction(globio_ssp3);
var ssp5_globio_constantUrban = globioUrbanMergeFunction(globio_ssp5);

//noRegain
var ssp1_globio_noRegain_constantUrban = globioUrbanMergeFunction(globio_ssp1_noRegain);
var ssp3_globio_noRegain_constantUrban = globioUrbanMergeFunction(globio_ssp3_noRegain);
var ssp5_globio_noRegain_constantUrban = globioUrbanMergeFunction(globio_ssp5_noRegain);


/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
//HABITAT SUITABLE RANGE ANALYSIS
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


// A function that calculates HSR change for each species
// Urban HSR change is calculated using the dedicated urban land forecasts
// Non-urban HSR change is calculated using the GLOBIO LULC forecasts
// Returns an image collection, with an image for each species. 
// Pixel values indicate proportion of pixel that represents HSR for species
// Bands for HSR in 2015 and 2050
// Output image properties contain sums of key statistics
// Inputs are:
// - "speciesRangeIC": Species expert range map, with habitat preferences as properties
// - "landuse": Reclassified GLOBIO LULC forecasts with urban land held constant from 2015 - 2050. Bands for 2015 and 2050.
// - "urbanProj": urban land projections for 2050
// - "globio_constantUrban": Image with urban land held constant from 2015 - 2050, but pasture and secondary vegetation are unchnged from the original GLOBIO LULC forescast.
// - "clusters": A feature collection of urban clusters 

var habitatProjectionFunction = function(speciesRangeIC, landUse, urbanProj, globio_constantUrban, clusters) {
  var hsrIC = speciesRangeIC.map(function(species) { //For each species
    //Get the species IGBP habitat preferences as a list
    var habitat_list = ee.List(species.get('habitats_lst'));
    //Create a list of the corresponding GLOBIO LULC categories
    //Details of habtat preference matching in Simkin et al. (2022) - Supplementary materials
    var remapList = habitat_list.map(function(num) {
      var remapDict = ee.Dictionary([
      '1', [70,71],
      '2', [50],
      '3', [80, 81, 82],
      '4', [60, 61],
      '5', [90],
      '6', [120],
      '7', [100, 110, 120, 121, 122, 130, 150],
      '8', [70, 80, 71, 60, 100, 72, 121],
      '9', [80, 62, 120, 70, 60, 100, 110],
      '10', [110, 120, 130, 140],
      '11', [160, 170, 180],
      '12', [2, 230, 231, 3, 4],
      '13', [1,190],
      '14', [2, 230, 231],
      '15', [220],
      '16', [151, 152, 153, 200, 201, 202]
      ]); //A dictionary of the prefs habitat category (keys) and corresponding GLOBIO LULC map categories (values)
      var remapDictKeys = ee.List(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16']); //Create a list of the keys
      var matchKey = remapDictKeys.get(ee.Number(num).subtract(1)); //For each number in 'habitat_list' return the corresponding key from 'remapDict'
      var match = remapDict.get(matchKey); //use the key to get the GLOBIO LULC category that matches the habitat preference
      return match; //return the list of habitat prefences
    }).flatten().distinct().sort();  //Final result is a list of matching GLOBIO LULC categories with duplicates removed in numeric order
    // Apply the reclassification
    // Returns a global image with all pixels that do not match habitat preferences masked
    // Pixel values are LULC categories
    var globio_habitatImg_2015 = landUse.remap({from: remapList, to: remapList, bandName: '2015'});
    var globio_habitatImg_2050 = landUse.remap({from: remapList, to: remapList, bandName: '2050'});

    // Merge into a single multi-band image
    var globio_habitatImg = ee.ImageCollection.fromImages([globio_habitatImg_2015, globio_habitatImg_2050]).toBands().select(['0_remapped', '1_remapped'], ['2015', '2050']);
    
    // Refine to the expert range map and elevation preferences of the species
    // Create an expert range map mask 
    // Species range: 1 = within range, 0 = outside range
    var region_mask = ee.Image(species).gt(0);
    
    // Create an elevation mask
    var elevation_mask = DEM.gte(ee.Number(species.get('elev_min'))).and(DEM.lte(ee.Number(species.get('elev_max'))));

    // Apply expert range map mask to the habitat image and add properties
    // Note that elevation mask is applied in later step
    var habitat = ee.Image(globio_habitatImg.updateMask(region_mask)
        .set('speciesName', species.get('name'))
        .set('taxa', species.get('taxa'))
        .set('prefs_habitat_list', habitat_list)
        .set('globio_habitat_list', remapList));
    
  
    // Create a geometry for species that do not have a geometry
    var projectionGeom = ee.Algorithms.If({
      condition: habitat.geometry().bounds().isUnbounded(),
      trueCase: ee.Geometry.Rectangle([-180,-85,180,85],'EPSG:4326',false),
      falseCase: habitat.geometry().bounds()
    });
    
    // Create a binary image habitat = 1, non-habitat = 0
    // note that the elevation mask has not yet been applied
    var habitatBinary = habitat.gt(0);//'2015' & '2050' bands


    //Create an image that shows proportion of HSR per pixel accounting for all forms of land use change
    // Non HSR pixels are masked
    // The value of HSR pixels not impacted by urban growth will be 1.
    // For species that do not have a habitat preference for urban land
    // - The value of HSR pixels impacted by urban growth will be equal to 1-the probability that the pixel that will become urban.
    //If species habitat preferences includes urban land:
    // - do not subtract the urban probability from the HSR pixels.
    // - add pixels to HSR image that become urban in 2050 if they are within the expert range map. The pixel value will be equal to the probability that the cell that becomes urban
    // Note that the elevation mask is applied here to avoid adding urban HSR back into areas that are outside of the elevation range of the species

    var habitatProps = ee.Image(ee.Algorithms.If({
        condition: ee.List(habitat.get('prefs_habitat_list')).contains(13),
        falseCase: habitatBinary.select('2015').addBands(habitatBinary.select('2050').subtract(urbanProj)) //From each pixel, subtract the proportion of the pixel predicted to become urban
                    .updateMask(elevation_mask), //Then apply elevation mask 
        trueCase: habitatBinary.select('2015').addBands(habitatBinary.select('2050').unmask(urbanProj).mask(species.mask())).selfMask() //If a pixel that is within the species range becomes urban, then add the urban value to the habitatBinary image, mask out any zeros
                                .updateMask(elevation_mask) //Then apply elevation mask
    }));//End of if statement
    
    // Area of HSR remaining after all forms of HSR loss (i.e. urban and non-urban)
    // Image with pixel value = pixel area
    var habitatArea = habitatProps.multiply(ee.Image.pixelArea());
    
    //Sum HSR area across species range
    //Gives total HSR area for 2015 and 2050 accounting for all forms of HSR loss
    var sumTable = habitatArea.reduceRegion({
        reducer: ee.Reducer.sum(),
        geometry: projectionGeom,
        scale: scale, 
        bestEffort:false,
        maxPixels:1e12
      });
    // Get the summed area of HSR for 2015 and 2050 as a number
    var hab2015 = ee.Number(sumTable.get('2015'));
    var hab2050 = ee.Number(sumTable.get('2050'));
 
    // An image that shows the change in habitat 2015 to 2050
    // Pixel values are equal to the proportion of habitat increase or decrease within the pixel
    // Pixel vaues: 0 = no change, +ve = increase, -ve = decrease
    // Zero values are masked so that only the change pixels are visible
    var changeImg_unmasked = habitatProps.select('2050').unmask(0).subtract(habitatProps.select('2015').unmask(0));
    var changeImg = changeImg_unmasked.updateMask(changeImg_unmasked.neq(0));
    
    // Find proportion of the change attributable to urban expansion
    // First, a layer that shows new urban land masked to the extent of the HSR pixels that show a change in HSR
    var newUrban = urbanProj.subtract(urban_2015).unmask(0).mask(changeImg.mask())  
    
    // Mask any pixels for which a species is predicted to both gain HSR due to non-urban land use change AND lose HSR due to urban land-use change (quite rare)
    // This change is counted as a net increase in HSR at the pixel level as part of the non-urban drivers section below, so mask them from this section
    // Create the mask
    var netGainUrbanLossMask = ee.Image(1).where(
      habitatBinary.select('2015').unmask(0).eq(0) //Note that habitatBinary has not been masked to elevation, so this may contain pixels outside the elevation range
        .and(habitatBinary.select('2050').eq(1))
        .and(habitatProps.select('2050').gt(0))//This step should remove any pixels outside the elevation range, as habitatProps is masked to the elevation range
        .and(habitatProps.select('2050').lt(1)),
      ee.Image(0) 
    );
    
    //Create a layer that shows the proportion of the change that is attributable to urban
    var urbanPropOfChange = ee.Image(1).subtract(changeImg.abs().subtract(newUrban))
      .updateMask(netGainUrbanLossMask) //mask out the net gain pixels as described above
    
    //Create an image that shows the area of HSR change attributable to urban
    //+ve areas indicate an increase, -ve values indicate a decrease
    var urbanChangeArea = ee.Image.pixelArea().multiply(changeImg)// area of all change
      .multiply(urbanPropOfChange).updateMask(netGainUrbanLossMask)
    
    //Get the summed area of HSR change due to urban land use change across the species range as a number  (meters)
    var urbanChange_m = urbanChangeArea.updateMask(urbanChangeArea.neq(0)).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: projectionGeom,
      scale: scale,
      bestEffort: false,
      maxPixels: 1e12
    }).getNumber('area');
    
    //Caluclate the remaining area of change - i.e. area of change attributable to non-urban drivers
    var nonUrbanLossArea = changeImg.multiply(ee.Image.pixelArea())// area of all change
      .multiply(ee.Image(1).subtract(urbanPropOfChange.unmask()))
    
    //Find the contribution of each of the non-urban drivers to HSR loss
    var reducerImg = nonUrbanLossArea.addBands(globio_constantUrban.select('2050')) //Adds the globio land use categories back as a band, mask out the zero values
      .updateMask(nonUrbanLossArea.neq(0)); //Mask out the pixels that have no change
    // Sum the area of HSR change caused by each LULC category
    var driversSum_m = reducerImg.reduceRegion({
        reducer: ee.Reducer.sum().group({
        groupField: 1,
        groupName: 'globioLC',
        }),
      geometry: projectionGeom,
      scale: scale,
      maxPixels: 1e12
    });
    
    // Get a list of the countries that contain HSR for the species in 2015 
    var presenceCountriesFC = countries.map(function(feat){
      var presenceAbsence = habitatProps.select('2015').reduceRegion({
        reducer: ee.Reducer.anyNonZero(),
        geometry: feat.geometry(),
        scale: scale,
        bestEffort: true
      }).getNumber('2015');
    return feat.set('presence', presenceAbsence)
    }).filter(ee.Filter.eq('presence', 1))
    
    //Calculate the area of HSR lost to urban in each country
    var urbanChangeArea_country = presenceCountriesFC.map(function(feat) {
      var changeArea = urbanChangeArea.reduceRegion({
        reducer: ee.Reducer.sum(),
        geometry:feat.geometry(),
        scale: scale,
        maxPixels: 1e12,
        bestEffort: false
      }).getNumber('area')
      return feat.set('urbanChangeArea_country', changeArea)
    }).toList(presenceCountriesFC.size().add(1))
      .map(function(i){return ee.Feature(i).toDictionary(['NAME', 'urbanChangeArea_country'])})
    
    //Add a list of the cities that contain HSR for the species in 2015
    var presenceCitiesFC = clusters.map(function(feat){
      var presenceAbsence = habitatProps.select('2015').reduceRegion({
        reducer: ee.Reducer.anyNonZero(),
        geometry: feat.geometry(),
        scale: scale,
        bestEffort: true
      }).getNumber('2015');
      var city = ee.Feature(feat).getString('cities');
      var cityIds = ee.Feature(feat).getString('cityIds');
    return feat.set('presence', presenceAbsence).set('city', city).set('cityIds', cityIds)
    }).filter(ee.Filter.eq('presence', 1))
    
    //Calculate net area of HSR change caused by each urban cluster (i.e. increase for species with habitat preference for urban or decrease otherwise)
    var urbanChangeArea_city = presenceCitiesFC.map(function(feat) {
      var changeArea = urbanChangeArea.reduceRegion({
        reducer: ee.Reducer.sum(),
        geometry:feat.geometry(),
        scale: scale,
        maxPixels: 1e12,
        bestEffort: false
      }).getNumber('area')
      return feat.set('urbanChangeArea_city', changeArea)
    }).toList(presenceCitiesFC.size().add(1))
      .map(function(feat) {return ee.Feature(feat).toDictionary(['city', 'urbanChangeArea_city', 'countryOfLargestCity', 'cityIds'])})
    
    //Get the area of urban change as a proportion of the 2015 HSR area, as a number
    //+ve = HSR gain, -ve = HSR loss
    var propUrbanChange =  urbanChange_m.divide(hab2015);
    
    //Get the area of all HSR change as a proportion of the 2015 HSR area, as a number
    var propAllChange = (hab2050.subtract(hab2015)).divide(hab2015)

    //species properties as a dictionary
    var properties = habitat.toDictionary();
    
    //Create an ee.Image to return
    //Multi-band image shows the HSR for the species in 2015 and 2050
    //pixel values indicate proportion of pixel that is suitible habitat for the species i.e. betwen 0 - 1
    //Add properties
    var Image = ee.Image(habitatProps).set(properties)
      .set('hab2015', hab2015)
      .set('hab2050', hab2050)
      .set('propUrbanChange', propUrbanChange)
      .set('propAllChange', propAllChange)
      .set('drivers_m', driversSum_m)
      .set('urbanChange_m', urbanChange_m)
      .set('urbanChangeArea_country', urbanChangeArea_country)
      .set('urbanChangeArea_city', urbanChangeArea_city);
      
    
    return Image;
    ////////////////////////////////////////////////////////////////////////////// 
    
  });
return hsrIC;
};


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//RESULTS AS AN IMAGE COLLECTION
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
var sampleSpp_ssp1_noRegain_results = ee.ImageCollection(habitatProjectionFunction(sampleSpp, ssp1_globio_noRegain_constantUrban, ssp1, ssp1_globio_constantUrban, ssp1Clusters)); 
print("sampleSpp_ssp1_noRegain_results", sampleSpp_ssp1_noRegain_results);
Map.setCenter(-91.25343155171775, 15.365242314016381, 7) 
Map.addLayer(sampleSpp_ssp1_noRegain_results.first().select(["2015"]), {min:0, max:1, palette:["red"]}, "2015_sampleSpp_ssp1_noRegain_results");
Map.addLayer(sampleSpp_ssp1_noRegain_results.first().select(["2050"]).selfMask(), {min:0, max:1, palette:["green"]}, "2050_sampleSpp_ssp1_noRegain_results");

var sampleSpp_ssp3_noRegain_results = ee.ImageCollection(habitatProjectionFunction(sampleSpp, ssp3_globio_noRegain_constantUrban, ssp3, ssp3_globio_constantUrban, ssp3Clusters)); 
var sampleSpp_ssp5_noRegain_results = ee.ImageCollection(habitatProjectionFunction(sampleSpp, ssp5_globio_noRegain_constantUrban, ssp5, ssp5_globio_constantUrban, ssp5Clusters));



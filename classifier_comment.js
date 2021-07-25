var randomPointsPreComputed = ee.FeatureCollection([geometry1,geometry2,geometry3,geometry4]);


//step0:numbers
var globOptions = { 
  bandSelect: ['green', 'swir1', 'swir2', 'nir', 'red'],
  bands2sr: ['B3', 'B11', 'B12', 'B8', 'B4'],
  bands2sr_TCI: ['TCI_G', 'B11', 'B12', 'B8', 'TCI_R'],
  bandaster: ['B01','B04','B09','B3N', 'B02'],
  startDate0: '2018-01-01',
  endDate0: '2018-12-31',
  startDate1: '2018-01-01',
  endDate1: '2018-03-31',
  startDate2: '2018-04-01',
  endDate2: '2018-06-30',
  startDate3: '2018-07-01',
  endDate3: '2018-09-30',
  startDate4: '2018-10-01',
  endDate4: '2018-12-31'
};





//step1:functions
var landsatFunctions = {
  //温度temperature
  //------------------------------
  // fromDN: function(image) {
  //   var bands = ['B10', 'B11', 'B12', 'B13', 'B14'] //thermal infrared
  //   var multiplier = ee.Image([0.006822, 0.006780, 0.006590, 0.005693, 0.005225])
  //   var k1 = ee.Image([3040.136402, 2482.375199, 1935.060183, 866.468575, 641.326517])
  //   var k2 = ee.Image([1735.337945, 1666.398761, 1585.420044, 1350.069147, 1271.221673])
  //   var radiance = image.select(bands).subtract(1).multiply(multiplier)
  //   var t = k2.divide(k1.divide(radiance).add(1).log())
  //   return t.select([0], ['temp']).copyProperties(image, ['system:time_start']);
  // },
  
  
  //------------------------------
  applySDWI: function(image){  
    var wi1 = image.expression(
    '(10*VV*VH)', //可能需要 abs(VV*VH)
    {
      'VV':image.select('VV'),
      'VH':image.select('VH')
    }
    )
                          .rename('Water1'); 
    var wi2 = wi1.log();
    var wi3 = wi2.expression(
    'water-8',
    {
      'water':wi2.select('Water1')
    }
    ).gt(0.2) //原本是0.35
                          .rename('sdwi');
    return image.addBands(wi3) .select('sdwi')
    
  },
  
  //------------------------------
  applyNDWI: function(image) {
    // apply NDWI to an image
    var ndwi = image.normalizedDifference(['green','nir']);
    return ndwi.select([0], ['ndwi']).copyProperties(image, ['system:time_start']);
  },
  
  //------------------------------
  applyMNDWI: function(image) {
  // apply MNDWI to an image
    var mndwi = image.normalizedDifference(["green","swir1"]);
    return image.select([0], ['mndwi']).copyProperties(image, ['system:time_start']);
  },
  
 //------------------------------
  applyNDVI: function(image) {
  // apply NDVI to an image
    var ndvi = image.normalizedDifference(['nir','red']);
    return ndvi.select([0], ['ndvi']).copyProperties(image, ['system:time_start']);
  },
 //------------------------------
  // applyAWEI: function(image) {
  //   // apply AWEI to an image
  //   var awei = image.expression("10*(b('green')-b('swir1'))-(0.25*b('nir')+2.75*b('swir2'))");
  //   return awei.select([0], ['awei']).copyProperties(image, ['system:time_start']);
  // },

  applyFMask: function maskS2clouds(image) {
     var qa = image.select('QA60');
     var cloudBitMask = ee.Number(2).pow(10).int();
     var cirrusBitMask = ee.Number(2).pow(11).int();
     var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and( qa.bitwiseAnd(cirrusBitMask).eq(0));
     return image.updateMask(mask).divide(10000);
  },

};









var reducer = ee.Reducer.min()
    .combine(ee.Reducer.max(), '', true)
    // .combine(ee.Reducer.stdDev(), '', true)
    // .combine(ee.Reducer.median(), '', true)
    //.combine(ee.Reducer.count(), '', true) 
    // .combine(ee.Reducer.percentile([10, 25, 50, 75,90]), '', true)
    // .combine(ee.Reducer.intervalMean(0, 10).setOutputs(['intMn0010']), '', true)
    // .combine(ee.Reducer.intervalMean(10, 25).setOutputs(['intMn1025']), '', true)
    // .combine(ee.Reducer.intervalMean(25, 50).setOutputs(['intMn2550']), '', true)
    // .combine(ee.Reducer.intervalMean(50, 75).setOutputs(['intMn5075']), '', true)
    // .combine(ee.Reducer.intervalMean(75, 90).setOutputs(['intMn7590']), '', true)
    // .combine(ee.Reducer.intervalMean(90, 100).setOutputs(['intMn90100']), '', true)
    // .combine(ee.Reducer.intervalMean(10, 90).setOutputs(['intMn1090']), '', true)
    // .combine(ee.Reducer.intervalMean(25, 75).setOutputs(['intMn2575']), '', true);


var geometry = ee.Geometry.Polygon(
        [[[121.729183, 31.664105],
          [121.729183, 31.264105],
          [122.129183, 31.264105],
          [122.129183, 31.664105]]], null, false);

// var randomPointsPreComputed = ee.FeatureCollection.randomPoints({
//     region: geometry,
//     points: 1000
//   });
  
  
  
  

var ASTER = ee.ImageCollection('ASTER/AST_L1T_003')
  .filterDate(globOptions.startDate0, globOptions.endDate0)
  .filterBounds(geometry);
  
// var ASTER_TEMP = ASTER.map(landsatFunctions.fromDN).aside(print);

var ASTER_1 = ASTER.select(globOptions.bandaster, globOptions.bandSelect)
.map(landsatFunctions.applyNDWI)
.reduce(reducer, globOptions.parallelScale)
.aside(print);  //NICE!

var ASTER_2 = ASTER.select(globOptions.bandaster, globOptions.bandSelect)
.map(landsatFunctions.applyMNDWI)
.reduce(reducer, globOptions.parallelScale)
.aside(print);  //NICE!

var ASTER_3 = ASTER.select(globOptions.bandaster, globOptions.bandSelect)
.map(landsatFunctions.applyNDVI)
.reduce(reducer, globOptions.parallelScale)
.aside(print);  //NICE!














var S1_GRD_VV1 = ee.ImageCollection("COPERNICUS/S1_GRD")
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
.filter(ee.Filter.eq('instrumentMode', 'IW'))
// .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
.select('VV')
.filterDate(globOptions.startDate1, globOptions.endDate1)
.filterBounds(geometry)
.map(function(image) {
  var edge = image.lt(-30.0);
  var maskedImage = image.mask().and(edge.not());
  return image.updateMask(maskedImage);
})
// .map(landsatFunctions.applySDWI)
// .reduce(reducer, globOptions.parallelScale)
.aside(print,'S1_GRD_VV1');

var S1_GRD_VV2 = ee.ImageCollection("COPERNICUS/S1_GRD")
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
.filter(ee.Filter.eq('instrumentMode', 'IW'))
.select('VV')
.filterDate(globOptions.startDate2, globOptions.endDate2)
.filterBounds(geometry)
.map(function(image) {
  var edge = image.lt(-30.0);
  var maskedImage = image.mask().and(edge.not());
  return image.updateMask(maskedImage);
})
// .map(landsatFunctions.applySDWI)
// .reduce(reducer, globOptions.parallelScale)
.aside(print,'S1_GRD_VV2');

var S1_GRD_VV3 = ee.ImageCollection("COPERNICUS/S1_GRD")
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
.filter(ee.Filter.eq('instrumentMode', 'IW'))
.select('VV')
.filterDate(globOptions.startDate3, globOptions.endDate3)
.filterBounds(geometry)
.map(function(image) {
  var edge = image.lt(-30.0);
  var maskedImage = image.mask().and(edge.not());
  return image.updateMask(maskedImage);
})
// .map(landsatFunctions.applySDWI)
// .reduce(reducer, globOptions.parallelScale)
.aside(print,'S1_GRD_VV3');

var S1_GRD_VV4 = ee.ImageCollection("COPERNICUS/S1_GRD")
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
.filter(ee.Filter.eq('instrumentMode', 'IW'))
.select('VV')
.filterDate(globOptions.startDate4, globOptions.endDate4)
.filterBounds(geometry)
.map(function(image) {
  var edge = image.lt(-30.0);
  var maskedImage = image.mask().and(edge.not());
  return image.updateMask(maskedImage);
})
// .map(landsatFunctions.applySDWI)
// .reduce(reducer, globOptions.parallelScale)
.aside(print,'S1_GRD_VV4');
  
var S1_GRD_VH1 = ee.ImageCollection("COPERNICUS/S1_GRD")
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
.filter(ee.Filter.eq('instrumentMode', 'IW'))
// .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
.select('VH')
.filterDate(globOptions.startDate1, globOptions.endDate1)
.filterBounds(geometry)
.map(function(image) {
  var edge = image.lt(-30.0);
  var maskedImage = image.mask().and(edge.not());
  return image.updateMask(maskedImage);
})
// .map(landsatFunctions.applySDWI)
// .reduce(reducer, globOptions.parallelScale)
.aside(print,'S1_GRD_VH1');

var S1_GRD_VH2 = ee.ImageCollection("COPERNICUS/S1_GRD")
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
.filter(ee.Filter.eq('instrumentMode', 'IW'))
.select('VH')
.filterDate(globOptions.startDate2, globOptions.endDate2)
.filterBounds(geometry)
.map(function(image) {
  var edge = image.lt(-30.0);
  var maskedImage = image.mask().and(edge.not());
  return image.updateMask(maskedImage);
})
// .map(landsatFunctions.applySDWI)
// .reduce(reducer, globOptions.parallelScale)
.aside(print,'S1_GRD_VH2');

var S1_GRD_VH3 = ee.ImageCollection("COPERNICUS/S1_GRD")
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
.filter(ee.Filter.eq('instrumentMode', 'IW'))
.select('VH')
.filterDate(globOptions.startDate3, globOptions.endDate3)
.filterBounds(geometry)
.map(function(image) {
  var edge = image.lt(-30.0);
  var maskedImage = image.mask().and(edge.not());
  return image.updateMask(maskedImage);
})
// .map(landsatFunctions.applySDWI)
// .reduce(reducer, globOptions.parallelScale)
.aside(print,'S1_GRD_VH3');

var S1_GRD_VH4 = ee.ImageCollection("COPERNICUS/S1_GRD")
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
.filter(ee.Filter.eq('instrumentMode', 'IW'))
.select('VH')
.filterDate(globOptions.startDate4, globOptions.endDate4)
.filterBounds(geometry)
.map(function(image) {
  var edge = image.lt(-30.0);
  var maskedImage = image.mask().and(edge.not());
  return image.updateMask(maskedImage);
})
// .map(landsatFunctions.applySDWI)
// .reduce(reducer, globOptions.parallelScale)
.aside(print,'S1_GRD_VH4');

//----------------------segmentation-----------------------------------------

var VV1 = ee.Image(S1_GRD_VV1.median()).clip(geometry);
var VV2 = ee.Image(S1_GRD_VV2.median()).clip(geometry);
var VV3 = ee.Image(S1_GRD_VV3.median()).clip(geometry);
var VV4 = ee.Image(S1_GRD_VV4.median()).clip(geometry);
var VH1 = ee.Image(S1_GRD_VH1.median()).clip(geometry);
var VH2 = ee.Image(S1_GRD_VH2.median()).clip(geometry);
var VH3 = ee.Image(S1_GRD_VH3.median()).clip(geometry);
var VH4 = ee.Image(S1_GRD_VH4.median()).clip(geometry);


var img_seg = VV1.addBands(VV2).addBands(VV3).addBands(VV4)
        .addBands(VV1).addBands(VV2).addBands(VV3).addBands(VV4).aside(print,'RGBVH+RGBVV')

var seeds = ee.Algorithms.Image.Segmentation.seedGrid(20).aside(print,'seeds');

var snic_pred = ee.Algorithms.Image.Segmentation.SNIC({
  image: img_seg, 
  compactness: 4,
  connectivity: 8,
  neighborhoodSize: 128,
  size: 10,
  seeds: seeds
}).aside(print,'SNIC')
    
var clusters_snic_pred = snic_pred.select("clusters")


var vectors_pred = clusters_snic_pred.reduceToVectors({
  geometryType: 'polygon',
  reducer: ee.Reducer.countEvery(),
  scale: 20,
  maxPixels: 1e13,
  geometry: geometry,
}); 

// var empty_pred = ee.Image().byte().aside(print,'byte');

// var outline_pred = empty_pred.paint({
//   featureCollection: vectors_pred,
//   color: 1,
//   width: 1
// }).aside(print,'paint');

// Map.addLayer(outline_pred, {palette: 'FF0000'}, 'vec_snic_pred');


//---------------- Classification ----------------------------

var FullImage = img_seg.toFloat().aside(print,'toFloat') ;

var train_points = randomPointsPreComputed.aside(print,'pretrainedpoints');

var train_polys = vectors_pred
.map(function(feat){
  feat = ee.Feature(feat);
  var point = feat.geometry();
//whether intersects
  var mappedPolys = train_points.map(function(poly){
    var cls = poly.get("Class")
    var intersects = poly.intersects(point, ee.ErrorMargin(1));
    var property = ee.String(ee.Algorithms.If(intersects, 'TRUE', 'FALSE'));
    return feat.set('belongsTo',  property).set('Class', cls);
  });
  return mappedPolys;
})
.flatten()
.filter(ee.Filter.neq('belongsTo', 'FALSE')).aside(print,'belongsToFALSEintersects');



var train_areas = train_polys //FeatureCollection (162 elements
  .reduceToImage({
    properties: ['Class'],
    reducer: ee.Reducer.first()
}).aside(print,'reduceToImage') //"first", double, EPSG:4326 Image (1 band)

// bands: List (1 element)
// 0: "first", double, EPSG:4326

.rename('Class').toInt().aside(print,'reduceToImagetoInt'); //"Class", signed int32, EPSG:4326 Image (1 band)

// bands: List (1 element)
// 0: "Class", signed int32, EPSG:4326


var vis_RF = {min: 0, max: 1,
palette: [ 'yellow' //0
,'green']//1
}

Map.addLayer(train_areas,vis_RF,"clip_img");


var predict_image = vectors_pred
  .reduceToImage({
    properties: ['label'],
    reducer: ee.Reducer.first()
}).rename('id').toInt();



FullImage = FullImage.addBands(predict_image).aside(print,'FullImage')//Image (7 bands)"id", signed int32, EPSG:4326

var FullImage_mean = FullImage.reduceConnectedComponents({
  reducer: ee.Reducer.mean(),
  labelBand: 'id'
});

var FullImage_max = FullImage.reduceConnectedComponents({
  reducer: ee.Reducer.max(),
  labelBand: 'id'
});

var FullImage_min = FullImage.reduceConnectedComponents({
  reducer: ee.Reducer.min(),
  labelBand: 'id'
});

var FullImage_std = FullImage.reduceConnectedComponents({
  reducer: ee.Reducer.stdDev(),
  labelBand: 'id'
});

var FullImage_median = FullImage.reduceConnectedComponents({
  reducer: ee.Reducer.median(),
  labelBand: 'id'
});

var Pred_bands = ee.Image.cat([
  FullImage_mean,
  FullImage_max,
  FullImage_std,
  FullImage_median,
  FullImage_min
]).float();

var predictionBands = Pred_bands.bandNames();

var training = Pred_bands.sample(train_points, 20).aside(print,'sample'); //FeatureCollection (159 elements, 34 columns)





var withRandom = training.randomColumn('random').aside(print,'randomColumn'); //FeatureCollection (159 elements, 35 columns)、




var split = 0.7;  // Roughly 70% training, 30% testing.
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));



var trainedClassifier = ee.Classifier.smileRandomForest({numberOfTrees:100}).train({
  features: trainingPartition,
  classProperty: 'Class',
  inputProperties: predictionBands
}).aside(print,'smileRandomForest'); //Classifier.train

var test = testingPartition.classify(trainedClassifier).aside(print,'classifytestingPartition'); //FeatureCollection (40 elements, 0 columns)


var classified_RF = Pred_bands.select(predictionBands).classify(trainedClassifier).aside(print,'classifyPred_bands'); //Image (1 band)  "classification", signed int32

var vis_RF = {min: 0, max: 1,
palette: [ 'yellow' //0
,'green']//1
}

Map.addLayer(classified_RF, vis_RF, 'Classified');




















 var S2_SR = ee.ImageCollection('COPERNICUS/S2_SR')
  .filterDate(globOptions.startDate0, globOptions.endDate0)
  .filterBounds(geometry)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 20)
  // .map(landsatFunctions.applyFMask)
  .select(globOptions.bands2sr_TCI, globOptions.bandSelect)
  // .sort('DATE_ACQUIRED',true)
  .aside(print,'S2_SR');
  
var S2_SR1 = S2_SR.map(landsatFunctions.applyNDWI)
.reduce(reducer, globOptions.parallelScale)
.aside(print); //NICE!

var S2_SR2 = S2_SR.map(landsatFunctions.applyMNDWI)
.reduce(reducer, globOptions.parallelScale)
.aside(print); //NICE!

var S2_SR3 = S2_SR.map(landsatFunctions.applyNDVI)
.reduce(reducer, globOptions.parallelScale)
.aside(print);  //NICE!

Map.setCenter(121.895238, 31.448851, 16); //崇明
// //显示sdwi
Map.addLayer(S2_SR1.select('ndwi_min'));  //NICE!
// //显示wendu
// var temperaturePalette = ['0571b0', '92c5de', 'f7f7f7', 'f4a582', 'ca0020']
// Map.addLayer(ASTER_1.select('temp'), {
//                   min: 273.15 - 20, max: 273.15 + 45, 
//                   palette: temperaturePalette
//                 }, 'ASTER temperature')

//STEP3:classify
// var trainComposite = S2_SR1
//       .addBands(S2_SR2)
//       .addBands(S2_SR3)
//       .addBands(S1_GRD)
//       .addBands(ASTER_1)
//       .addBands(ASTER_2)
//       .addBands(ASTER_3).randomColumn('random', 2015).aside(print); 
// var bands = trainComposite.bandNames();


// var train = trainComposite.randomColumn('random', 2015); //seed = 2015

// var trainingSet = trainComposite.filter(ee.Filter.gte('random',globOptions.trainingValidationRatio)).aside(print);
// var validationSet = trainComposite.filter(ee.Filter.lt('random',globOptions.trainingValidationRatio)).aside(print);


//4. Training Data
// var predictorSet = randomPointsPreComputed.aside(print);

// 5. Classify 

// ee.Classifier
//   classProperty: 'landcover',
//   inputProperties: bands
// });

// var classifier = ee.Classifier.smileRandomForest(globOptions.nTrees)
//   .train({
//     features: trainingSet,
//     classProperty: 'CLASS',
//     inputProperties: bands
//   })
//   .setOutputMode('CLASSIFICATION').aside(print);
  
// var classified = trainComposite.select(bands)
//   .classify(classifier).aside(print);
  
// // 海岸区域界定
// var  finalMask =  ee.Image('UQ/murray/Intertidal/v1_1/data_mask');
// var finalOut = classified.byte()
//   .mask(finalMask).aside(print); 

var randomPointsPreComputed = ee.FeatureCollection([geometry1,geometry2,geometry3,geometry4]);


//step0:numbers
var globOptions = { 
  bandSelect: ['green', 'swir1', 'swir2', 'nir', 'red'],
  bands2sr: ['B3', 'B11', 'B12', 'B8', 'B4'],
  bands2sr_TCI: ['TCI_G', 'B11', 'B12', 'B8', 'TCI_R'],
  bandaster: ['B01','B04','B09','B3N', 'B02'],
  startDate1: '2016-01-01',
  endDate1: '2016-03-31',
  startDate2: '2016-04-01',
  endDate2: '2016-06-31',
  startDate3: '2016-07-01',
  endDate3: '2016-09-31',
  startDate4: '2016-10-01',
  endDate4: '2016-12-31'
};




//step1:functions
var landsatFunctions = {
  //温度temperature
  //------------------------------
  fromDN: function(image) {
    var bands = ['B10', 'B11', 'B12', 'B13', 'B14'] //thermal infrared
    var multiplier = ee.Image([0.006822, 0.006780, 0.006590, 0.005693, 0.005225])
    var k1 = ee.Image([3040.136402, 2482.375199, 1935.060183, 866.468575, 641.326517])
    var k2 = ee.Image([1735.337945, 1666.398761, 1585.420044, 1350.069147, 1271.221673])
    var radiance = image.select(bands).subtract(1).multiply(multiplier)
    var t = k2.divide(k1.divide(radiance).add(1).log())
    return t.select([0], ['temp']).copyProperties(image, ['system:time_start']);
  },
  
  
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

var randomPointsPreComputed = ee.FeatureCollection.randomPoints({
    region: geometry,
    points: 1000
  });

var ASTER = ee.ImageCollection('ASTER/AST_L1T_003')
  .filterDate(globOptions.startDate, globOptions.endDate)
  .filterBounds(geometry)
  
var ASTER_TEMP = ASTER.map(landsatFunctions.fromDN).aside(print);

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


var S1_GRD = ee.ImageCollection("COPERNICUS/S1_GRD")
//             .filter(ee.Filter.eq('instrumentMode', 'IW'))
.filterDate(globOptions.startDate, globOptions.endDate)
.filterBounds(geometry)
.map(function(image) {
  var edge = image.lt(-30.0);
  var maskedImage = image.mask().and(edge.not());
  return image.updateMask(maskedImage);
})
.map(landsatFunctions.applySDWI)
.reduce(reducer, globOptions.parallelScale)
.aside(print);
  

 var S2_SR = ee.ImageCollection('COPERNICUS/S2_SR')
  .filterDate(globOptions.startDate, globOptions.endDate)
  .filterBounds(geometry)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 20)
  .map(landsatFunctions.applyFMask)
  .select(globOptions.bands2sr_TCI, globOptions.bandSelect)
  .sort('DATE_ACQUIRED',true).aside(print);
  
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
var trainComposite = S2_SR1
      .addBands(S2_SR2)
      .addBands(S2_SR3)
      .addBands(S1_GRD)
      .addBands(ASTER_1)
      .addBands(ASTER_2)
      .addBands(ASTER_3).randomColumn('random', 2015).aside(print); 
var bands = trainComposite.bandNames();


// var train = trainComposite.randomColumn('random', 2015); //seed = 2015

var trainingSet = trainComposite.filter(ee.Filter.gte('random',globOptions.trainingValidationRatio)).aside(print);
var validationSet = trainComposite.filter(ee.Filter.lt('random',globOptions.trainingValidationRatio)).aside(print);


//4. Training Data
var predictorSet = randomPointsPreComputed.aside(print);

// 5. Classify 

// ee.Classifier
//   classProperty: 'landcover',
//   inputProperties: bands
// });

var classifier = ee.Classifier.smileRandomForest(globOptions.nTrees)
  .train({
    features: trainingSet,
    classProperty: 'CLASS',
    inputProperties: bands
  })
  .setOutputMode('CLASSIFICATION').aside(print);
  
var classified = trainComposite.select(bands)
  .classify(classifier).aside(print);
  
// 海岸区域界定
var  finalMask =  ee.Image('UQ/murray/Intertidal/v1_1/data_mask');
var finalOut = classified.byte()
  .mask(finalMask).aside(print); 

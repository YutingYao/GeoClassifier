

//////////////////////////////////////////////////////////////////////////////////////// 
// geometry1 - geometry4 是标注好的数据集，采用监督分类的方法，所以需要预先对数据进行手动分类
// table1 - table44 是海岸线，把shp格式的海岸线边界导入GEE
// Developed in Google Earth Engine:
// https://code.earthengine.google.com
//////////////////////////////////////////////////////////////////////////////////////// 
var site = ee.Geometry.Polygon([-180, 60, 0, 60, 180, 60, 180, -60, 10, -60, -180, -60], null, false);


var globOptions = { 
  bandSelect: ['green',  'nir', 'red'],
  startDate0: '2018-01-01',
  endDate0: '2020-12-31',
  startDate1: '2018-01-01',
  endDate1: '2018-03-31',
  startDate2: '2018-04-01',
  endDate2: '2018-06-30',
  startDate3: '2018-07-01',
  endDate3: '2018-09-30',
  startDate4: '2018-10-01',
  endDate4: '2018-12-31',
};

var season1 = ee.Filter.date(globOptions.startDate1, globOptions.endDate1);
var season2 = ee.Filter.date(globOptions.startDate2, globOptions.endDate2);
var season3 = ee.Filter.date(globOptions.startDate3, globOptions.endDate3);
var season4 = ee.Filter.date(globOptions.startDate4, globOptions.endDate4);

var tablemerge = table1.merge(table2).merge(table3).merge(table4).merge(table5);
var tablemerge_af = table6.merge(table7);
var tablemerge_an = table8.merge(table9);
var tablemerge_eu = table10.merge(table11).merge(table12).merge(table13).merge(table14).merge(table15).merge(table16).merge(table17);
var tablemerge_na = table18.merge(table19).merge(table20).merge(table21).merge(table22).merge(table23).merge(table24).merge(table25);
                       .merge(table26).merge(table27).merge(table28).merge(table29).merge(table30).merge(table31).merge(table32).merge(table33).merge(table34).merge(table35);
var tablemerge_oa = table36.merge(table37).merge(table38);
var tablemerge_sa = table39.merge(table40).merge(table41).merge(table42).merge(table43).merge(table44);

var randomPointsPreComputed = ee.FeatureCollection([geometry1,geometry2,geometry3,geometry4]);

// 1. Functions

var wetlandFunctions = {

  applyFMask: function maskS2clouds(image) {
    var qa = image.select('QA60');
    var cloudBitMask = ee.Number(2).pow(10).int();
    var cirrusBitMask = ee.Number(2).pow(11).int();
    var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and( qa.bitwiseAnd(cirrusBitMask).eq(0));
    return image.updateMask(mask).divide(10000);
  },

  applyNDWI: function(image) {
    // apply NDWI to an image
    var ndwi = image.normalizedDifference(['green','nir']);
    return ndwi.select([0], ['ndwi']);
  },



  applyNDVI: function(image) {
    // apply NDVI to an image
    var ndvi = image.normalizedDifference(['nir','red']);
    return ndvi.select([0], ['ndvi']);
  },
  
  applyNeigh: function(image){
    return image.reduceNeighborhood({
  reducer: ee.Reducer.mean(),
  kernel: ee.Kernel.circle(4),
  });
  },
  
  applySDWI: function(image){
    var wi1 =image.expression(
      '(10*abs(VV*VH))', //可能需要 abs(S1mean*VH)
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
      ).gt(0.35) //原本是0.35
                            .rename('sdwi');
    return image.addBands(wi3).select('sdwi');
  }


  
};

///////////////////////////////////////////
// 计算visual
/////////////////////////////////////////
var dataset = ee.ImageCollection('COPERNICUS/S2_SR')
                  .filterDate('2021-01-01', '2021-01-30')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20))
                  .map(wetlandFunctions.applyFMask);

var visualization = {
  min: 0.0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2'],
};

Map.addLayer(dataset.mean(), visualization, 'RGB');
///////////////////////////////////////////
// 计算visual
/////////////////////////////////////////

var S1_1 = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(season1)
        .filter(ee.Filter.intersects('.geo', tablemerge.geometry(), null, null, 1000))
        .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
        .filterMetadata('resolution_meters', 'equals' , 10)
        .select(['VV',  'VH']);

var S1_2 = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(season2)
        .filter(ee.Filter.intersects('.geo', tablemerge.geometry(), null, null, 1000))
        .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
        .filterMetadata('resolution_meters', 'equals' , 10)
        .select(['VV',  'VH']);

var S1_3 = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(season3)
        .filter(ee.Filter.intersects('.geo', tablemerge.geometry(), null, null, 1000))
        .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
        .filterMetadata('resolution_meters', 'equals' , 10)
        .select(['VV',  'VH']);
                
var S1_4 = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(season4)
        .filter(ee.Filter.intersects('.geo', tablemerge.geometry(), null, null, 1000))
        .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
        .filterMetadata('resolution_meters', 'equals' , 10)
        .select(['VV',  'VH']);

print(S1_1.size(),S1_2.size(),S1_3.size(),S1_4.size());


var S2collection = ee.ImageCollection('COPERNICUS/S2_SR')
                 .filterDate(globOptions.startDate0, globOptions.endDate0)
                 .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 10)
                 .filter(ee.Filter.intersects('.geo', tablemerge.geometry(), null, null, 1000))
                                  .map(wetlandFunctions.applyFMask)
                 .select(['B3',  'B8', 'B4'], globOptions.bandSelect);

var S2collectionb = ee.ImageCollection('COPERNICUS/S2_SR')
                 .filterDate(globOptions.startDate0, globOptions.endDate0)
                 .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 10)
                 .filter(ee.Filter.intersects('.geo', tablemerge.geometry(), null, null, 1000))
                                  .map(wetlandFunctions.applyFMask)
                 .select('B4','B5','B12');

var swirmean = S2collectionb.map(wetlandFunctions.applyNeigh).mean();
// 0.9581147763919647

var ndwimax = S2collection.map(wetlandFunctions.applyNDWI).map(wetlandFunctions.applyNeigh).max();
var ndwimin = S2collection.map(wetlandFunctions.applyNDWI).map(wetlandFunctions.applyNeigh).min();
var ndwistd = S2collection.map(wetlandFunctions.applyNDWI).map(wetlandFunctions.applyNeigh).reduce(ee.Reducer.stdDev());

                            
var ndvimax = S2collection.map(wetlandFunctions.applyNDVI).map(wetlandFunctions.applyNeigh).max();
var ndvimin = S2collection.map(wetlandFunctions.applyNDVI).map(wetlandFunctions.applyNeigh).min();
var ndvistd = S2collection.map(wetlandFunctions.applyNDVI).map(wetlandFunctions.applyNeigh).reduce(ee.Reducer.stdDev());


var s1max = S1_1.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).max();
var s1min = S1_1.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).min();
var s1std = S1_1.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).reduce(ee.Reducer.stdDev());//.reduce(ee.Reducer.stdDev())


var s2max = S1_2.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).max();
var s2min = S1_2.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).min();
var s2std = S1_2.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).reduce(ee.Reducer.stdDev());


var s3max = S1_3.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).max();
var s3min = S1_3.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).min();
var s3std = S1_3.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).reduce(ee.Reducer.stdDev());

var s4max = S1_4.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).max();
var s4min = S1_4.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).min();
var s4std = S1_4.map(wetlandFunctions.applySDWI).map(wetlandFunctions.applyNeigh).reduce(ee.Reducer.stdDev());



var FullImage = ee.Image.cat(swirmean,ndwimax,ndwimin,ndwistd,ndvimax,ndvimin,ndvistd,s1max,s1min,s1std,s2max,s2min,s2std,s3max,s3min,s3std,s4max,s4min,s4std).aside(print); //


var bandNames = FullImage.bandNames();



// 4. Training Data

var train_points = randomPointsPreComputed;

var datasamples = FullImage.select(bandNames).sampleRegions({
  collection: train_points,
  properties: ['Class'],
  scale: 10
});


var withRandom = datasamples.randomColumn('random').aside(print,'randomColumn'); //FeatureCollection (159 elements, 35 columns)、
var split = 0.7;  // Roughly 70% training, 30% testing.
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));



var trainedClassifier = ee.Classifier.smileRandomForest({numberOfTrees:100}).train({
  features: trainingPartition,
  classProperty: 'Class',
  inputProperties: bandNames
}).aside(print,'smileRandomForest'); //Classifier.train

var test = testingPartition.classify(trainedClassifier)
        .aside(print,'classifytestingPartition'); //FeatureCollection (40 elements, 0 columns)


var classified_RF = FullImage.select(bandNames).classify(trainedClassifier)
        .aside(print,'classifyPred_bands'); //Image (1 band)  "classification", signed int32

var vis_RF = {min: 0, max: 3, palette: ['blue','Cyan','green', 'red']};

Map.addLayer(classified_RF, vis_RF, 'Classified',false);


// Extra post-process

//Run the classification for test data
var validation = testingPartition.classify(trainedClassifier);
var errorMatrix = validation.errorMatrix('Class', 'classification');
print('SENTINEL Error Matrix:', errorMatrix);
print('SENTINEL SVM Accuracy:', errorMatrix.accuracy());  
print('SENTINEL SVM Comsumers Accuracy:', errorMatrix.consumersAccuracy());
print('SENTINEL SVM Producers Accuracy:', errorMatrix.producersAccuracy());
print('SENTINEL SVM Kappa:', errorMatrix.kappa());

var mud = classified_RF.eq(1).selfMask().aside(print);

Export.image.toAsset({
  image: mud, 
  description: 'mud30',
  region: site, 
  maxPixels: 1e13,
  crs: "EPSG:4326",
  scale: 30,

});





Map.addLayer(mud,{palette: ['yellow']},'mud');

var veg = classified_RF.eq(2).selfMask().aside(print);

Map.addLayer(veg,{palette: ['green']},'tree',false);

var other = classified_RF.eq(3).selfMask().aside(print);

Map.addLayer(other,{palette: ['red']},'others');


// ------------------------------------------------------------------
Map.setCenter(120.8593, 32.0474, 8);  //南通
Map.setCenter(119.4382, 37.1148, 12); //SHANDONG
Map.setCenter(121.262, 30.3646, 12); //qiantangjiang

Map.setCenter(121.3894, 30.4808, 12); //qiantangjiang
Map.setCenter(139.0549328, -6.8820491, 12); //mappi regency
Map.setCenter(110.244291, 20.032579, 17); // , 西秀
Map.setCenter(121.895238, 31.448851, 14); //崇明
Map.setCenter(121.97263, 31.47507, 14); //崇明2


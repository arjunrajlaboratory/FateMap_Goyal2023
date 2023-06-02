%Test data can be found at https://www.dropbox.com/sh/lrg0v3y6o28on4q/AABB1Z2NXTVUsaKO0NXJBnepa?dl=0
test = scanObject('scanFile', '20190721_100908_528__ChannelCY3,A594,CY5,DAPI_Seq0000.nd2', 'tilesTable', 'tilesTable.csv');

scanDim = [50 50];
scanMatrix = vec2mat(1:scanDim(1)*scanDim(2), scanDim(1));
for i = 2:2:scanDim(1)
    scanMatrix(i, :) = fliplr(scanMatrix(i, :));
end

test.scanMatrix = scanMatrix;
test.stitchDAPI(scanMatrix);
test.contrastDAPIstitch();
test.loadStitches();
test.resizeStitchedScans();

testMask = maskTable(test, 'masksTest.csv');

testNuclei = nucleiTable(test, testMask);
testNuclei.stitchDAPImask();
testNuclei.findNuclei();
testNuclei.addColors();
testSpots = spotTable(test, testMask, testNuclei, 'spotsFile', 'spots.csv');
testSpots.assignSpotsToNuclei();

testSpots.makeCentroidList();
testSpots.defaultThresholds();

%% To launch GUI 
p = d2ThresholdView2(test, testSpots, testMask);
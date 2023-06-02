
% cd('/Users/arjunraj/Downloads/iCards_positive/w2');

ft = fileTable;
ft.loadFilesLegacy();

ct = cellTable();
ct.findCells(ft);

st = spotTable(ft);
st.findSpots(ft);
ints = st.getIntensities("tmr");
histogram(ints);
st.applyThreshold("tmr",17);


st.assignSpotsToCells(ct);

rect = [50 50 1400 1400];
outSpots = st.getSpotsInRect("tmr",rect);
[im,outRect] = ft.getImageFromRect(rect,"tmr",'largest');
imshow(im,[])
hold on;
[x,y] = d2utils.globalToLocalCoords(outRect,outSpots.x,outSpots.y);
%scatter(y,x,[],outSpots.nearestCellID)
plot(y,x,'g*');


figure
colormap prism
scatter(y,x,[],outSpots.nearestCellID)

currCells = ct.getCellsInRect(rect);
[x2,y2] = d2utils.globalToLocalCoords(outRect,currCells.x,currCells.y);
hold on;
plot(y2,x2,'c*');

% To mask points that are in h:
h = drawfreehand;
st.addMask(h.Position,outRect,"tmr");
% To get all non-masked points.
st.maskAllPoints();
outSpots = st.getValidNonmaskedSpotsInRect("tmr",rect);


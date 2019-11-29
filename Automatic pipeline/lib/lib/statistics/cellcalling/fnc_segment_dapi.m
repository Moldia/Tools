function [CellMap, DapiBoundaries] = fnc_segment_dapi(cropped_cell_image, threshhold)
% [CellMap, DapiBoundaries] = fnc_segment_dapi(cropped_cell_image, threshhold)
% second input optional
% segment dapi, adapted from iss suite by Kenneth Harris
% Xiaoyan, 2018

%% adjust and threshold
Dapi = imadjust(cropped_cell_image); % contrast enhancement
ImSz = size(Dapi);

if nargin < 2
    ThreshVal = prctile(Dapi(:), 90);
else
    ThreshVal = threshold;
end

bwDapi = imerode(Dapi>ThreshVal, strel('disk', 2));

%% find local maxima 
dist = bwdist(~bwDapi);
dist0 = dist;
% remove objects smaller than 5 px in radius
dist0(dist<5)=0;
ddist = imdilate(dist0, strel('disk', 7));
%clear dist 
impim = imimposemin(-dist0, imregionalmax(ddist));
clear dist0

figure(301);
subplot(2,1,1)
imagesc(dist);
subplot(2,1,2)
imagesc(impim);

%% segment
% remove pixels at watershed boundaries
bwDapi0 = bwDapi;
bwDapi0(watershed(impim)==0)=0;

% assign all pixels a label
labels = uint32(bwlabel(bwDapi0));
[d, idx] = bwdist(bwDapi0);

% now expand the regions by a margin
CellMap0 = zeros(ImSz, 'uint32');
Expansions = (d<10);
CellMap0(Expansions) = labels(idx(Expansions));

% get rid of cells that are too small
rProps0 = regionprops(CellMap0); % returns XY coordinate and area
BigEnough = [rProps0.Area]>200;
NewNumber = zeros(length(rProps0),1);
NewNumber(~BigEnough) = 0;
NewNumber(BigEnough) = 1:sum(BigEnough);
CellMap = CellMap0;
CellMap(CellMap0>0) = NewNumber(CellMap0(CellMap0>0));

figure(302)
image(label2rgb(CellMap, 'jet', 'w', 'shuffle'));


%%
% CellYX = fliplr(vertcat(rProps(BigEnough).Centroid)); % because XY

%% make image with boundaries
Boundaries = (CellMap ~= imdilate(CellMap,strel('disk', 1)));
DapiBoundaries = Dapi;

OddPix = mod((1:size(Dapi,2)) + (1:size(Dapi,1))', 2);
DapiBoundaries(Boundaries & OddPix) = .3 * max(Dapi(:));
DapiBoundaries(Boundaries & ~OddPix) = 0;

end


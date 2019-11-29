% Ab detection
% test pipeline for one image set (tile or full image)
% Xiaoyan, 2017

% reader = bfGetReader('F:\161220KI\161230_161220KI_3-1_Ab.czi');
% reader.setSeries(43-1);
% iPlane = reader.getIndex(5-1, 1-1, 0)+1;
% imDapi = bfGetPlane(reader, iPlane);
% iPlane = reader.getIndex(5-1, 3-1, 0)+1;
% imAb = bfGetPlane(reader, iPlane);

%% DAPI
imDapi = imread('E:\Melanoma\Align\Base_1_aligned-1.tif');

blurDapi = imfilter(imDapi, fspecial('average', 5));
% threshold
% bwDapi = im2bw(blurDapi, graythresh(blurDapi));  
bwDapi = im2bw(blurDapi, .01);   % manual threshold
bwDapi = imfill(bwDapi, 'holes');

% suppress local maxima and force seed positions (shape based)
distDapi = bwdist(~bwDapi);                 % distance transform
maxmask = getnhood(strel('disk', 7));       % maxima suppression size
maximaDapi = ones(size(distDapi));
maxfilter = ordfilt2(distDapi, sum(maxmask(:)), maxmask);   % maximum filter
maximaDapi(distDapi < maxfilter) = 0;
imposedDapi = imimposemin(-distDapi, maximaDapi);   % impose minima, force seed positions

% % suppress local maxima and force seed positions (intensity based)
% maxmask = getnhood(strel('disk', 15));    % maxima suppression size (big because of nucleoli)
% maximaDapi = blurDapi;
% maximaDapi(blurDapi < ordfilt2(blurDapi, sum(maxmask(:)), maxmask)) = 0;    % maximum filter
% maximaDapi(~bwDapi) = 0;    % remove outside objects
% imposedDapi = imimposemin(65535-blurDapi, maximaDapi);

% watershed to separate touching objects
wsBoundariesDapi = watershed(imposedDapi) > 0;
bwDapi = bwDapi & wsBoundariesDapi;

% size threshold (between 8 and 100 in diameter)
propDapi = struct2cell(regionprops(bwDapi, 'EquivDiameter', 'PixelIdxList'));
excludeDapi = cellfun(@(v) sum(maximaDapi(v))==0, propDapi(2,:));       % no maxima
excludeDapi = excludeDapi | ~cellfun(@(v) v>8 & v<100, propDapi(1,:));  % size threshold

% label image
labelDapi = bwlabel(bwDapi);
labelDapi(ismember(labelDapi, find(excludeDapi))) = 0;
bwDapi = logical(labelDapi);
[~, ~, labelDapi] = unique(labelDapi);  % relabel
labelDapi = reshape(labelDapi-1, size(imDapi));
outlinesDapi = bwDapi - imerode(bwDapi, strel('disk',1));
% imshow(repmat(imDapi*10,1,1,3) + cat(3, uint16(outlinesDapi)*65535, zeros(size(imDapi)), zeros(size(imDapi))))

clearvars -except bwDapi labelDapi outlinesDapi

%% Kv2.1
imAb = imread('E:\Melanoma\Align\Base_3_aligned-1.tif');

blurAb = imfilter(imAb, fspecial('average', 5));
bwAb = im2bw(blurAb, .01);
bwAb = imfill(bwAb, 'holes');
bwAb = imopen(bwAb, strel('disk', 7));

% Sobel filter
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(blurAb), hy, 'replicate');
Ix = imfilter(double(blurAb), hx, 'replicate');
sobelAb = sqrt(Ix.^2 + Iy.^2);

% watershed on Sobel filtered image
markerAb = bwDapi | ~bwAb;      % force Ab extends at least as far as DAPI
outlines = imdilate(bwDapi, ones(3)) - bwDapi;
markerAb(outlines==1) = 0;      % keep separated DAPI
shrunkAb = bwmorph(bwAb, 'shrink', inf);  % central pixel of all objects in bwAb, in order to keep all Ab objects
imposedAb = imimposemin(sobelAb, markerAb|shrunkAb);
wsAb = watershed(imposedAb);    

% label Ab object image
labelAb = logical(wsAb);        % wsAb: watershed lines=0, obejcts>=1
labelAb = bwlabel(labelAb);
labelAb = labelAb - bwareaopen(labelAb, 5000);  % remove super big objects, i.e. background
bwAb_new = logical(labelAb);

bwAb_new = imfill(bwAb_new, 'holes');
outlinesAb = bwAb_new - imerode(bwAb_new, strel('disk',1));
% imshow(repmat(imAb*20,1,1,3) + cat(3, uint16(outlinesDapi)*65535, uint16(outlinesAb)*65535, zeros(size(imDapi))))

% UNCOMMENT HERE IF WANT TO KEEP ONLY DAPI-POSITIVE CELLS
% % relabel according to DAPI objects, will also remove Ab objects not
% % overlapping with DAPI
% labelAb_new = bwlabel(bwAb_new);
% [DapiLabels, DapiLabelLocations] = unique(labelDapi);
% AbLabels = labelAb_new(DapiLabelLocations(2:end));  % skip label=0
% [~, labelAb_new] = ismember(labelAb_new, AbLabels);
% % temp = floor(rand(1)*max(labelAb_new(:)));
% % imshow(cat(3, labelDapi==temp, labelAb_new==temp, zeros(size(imAb))))

% UNCOMMENT HERE IF WANT TO KEEP ALL CELLS
% label Ab objects and find corresponding DAPI objects
labelAb_new = bwlabel(bwAb_new);
labelLocations = find(labelAb_new);
labels = labelAb_new(labelLocations);
mapAb_Dapi = [0 full(max(sparse(labelLocations, labels, labelDapi(labelLocations))))];
% labelAb_Dapi = mapAb_Dapi(labelAB_new+1);
% temp = floor(rand(1)*max(labelAb_new(:)));
% imshow(double(cat(3, labelDapi==temp, labelAb_Dapi==temp, labelAb_new==temp)))

% expand 10 pixels
[D, labelExpandedAb] = bwdist(bwAb_new);
labelExpandedAb = double(labelExpandedAb);
labelExpandedAb = labelAb_new(labelExpandedAb);
labelExpandedAb(D>10) = 0;
outlinesExpandedAb = labelExpandedAb ~= imerode(labelExpandedAb, strel('disk',1));
outlines = cat(3, uint16(outlinesDapi)*65535, uint16(outlinesExpandedAb)*65535, zeros(size(imDapi)));
imshow(repmat(imAb*20,1,1,3) + outlines)

% imwrite(outlines, 'outlines.png');

%% cell centroid
propCell = regionprops(labelExpandedAb, 'Centroid');    % or use labelAb_new instead
centroidCell = cat(1,propCell.Centroid);

% save before?
clearvars -except labelExpandedAb centroidCell

%% relate spots to cells (check pixel overlap, not precise)
% NOT TESTED
% get spot coordinates from fig file
scale = .25;
ch = get(gca, 'child');
total = 0;
for i = 1:length(ch)
    if strcmp(ch(i).Type, 'line')
        total = [total; total(end)+length(ch(i).XData)];
    end
end

pos = zeros(total(end),2);
nameidx = zeros(total(end),1);
names = [];
for i = 1:length(ch)-1
    pos(total(i)+1:total(i+1),:) = [ch(i).XData', ch(i).YData'];
    nameidx(total(i)+1:total(i+1)) = i;
    names = [names; {ch(i).DisplayName}];
end
pos = pos/scale;

% find parent cell
roundedCoord = round(pos);
linearCoord = (roundedCoord(:,1)-1)*size(labelExpandedAb,1) + roundedCoord(:,2);
idCell = labelExpandedAb(linearCoord);      % unlikely to have extreme values (e.g. spot rounded to coordinate(0,0)), should be fine

perCellCount = hist3([idCell, nameidx],...
    [{unique(labelExpandedAb)}, {1:length(names)}]);

% save

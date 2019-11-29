function [imLabel, objCentroid, objProps, objPixel] = ...
    filterobjects(imBW, threshProp, threshold)
% find connected components in a binary image
% filter objects based on threshProp property
% Xiaoyan, 2017


imCC = bwconncomp(imBW, 8);
imLabel = labelmatrix(imCC);
objProps = regionprops(imCC, 'Centroid', threshProp);
objCentroid = cat(1, objProps.Centroid);
objProps = cat(1, objProps.(threshProp));
objPixel = imCC.PixelIdxList;

% threshold
toExclude = find(objProps < threshold);
imLabel(ismember(imLabel, toExclude)) = 0;
objCentroid(toExclude,:) = [];
objProps(toExclude,:) = [];
objPixel(toExclude) = [];

% re-order
[~, ~, imLabel] = unique(imLabel);
imLabel = reshape(imLabel, size(imBW,1), size(imBW,2)) - 1;

end
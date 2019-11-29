function [grid_shiftright, grid_shiftleft, xypos] = hexbin_coord2grid(xypos, radius)
% [grid_shiftright, grid_shiftleft, xypos] = hexbin_coord2grid(xypos, radius)
% convert hexbin centroid coordinates back to grid
% returns two arrays, right and left shifted
% Xiaoyan, 2017

% expected distance between centroids
dy = 3/2 * radius;
dx = sqrt(3) * radius;


% number of grids in X and Y direction
nX = floor(round(range(xypos(:,1))/dx)) + 1;
nY = floor(round(range(xypos(:,2))/dy)) + 1;

assert(nX&nY, 'number of hexbins cannot be zero in any  direction');

% convert to grid coordinates
xypos = bsxfun(@minus, xypos, min(xypos, [], 1));
xypos = bsxfun(@rdivide, xypos, [dx, dy]) + 1;
xypos = round(xypos, 1);

% accumulated arrays
grid_shiftright = accumarray(fliplr(ceil(xypos)), 1, [nY nX+1]);
grid_shiftleft = accumarray(fliplr(floor(xypos)), 1, [nY nX+1]);

end

function [tilePosYX, gridpos] = create_tilepos_from_image(image, imscale,...
    imorder, tilesize, overlap)
% [tilepos, gridpos] = create_tilepos_from_image(image, imscale,...
%     imorder, tilesize, overlap)
% create tile position from a stitched image when metadata is corrupted
% Xiaoyan, 2017


if nargin < 2
    imscale = 1;
end

if nargin < 3
    % default: snake
    imorder = 'EWS';
end

if nargin < 4
    % default image size: 2048 x 2048 px
    tilesize = 2048;
end

if nargin < 5
    % default tile overlap in imaging: 10%
    overlap = .1;
end

unitsize = tilesize*(1-overlap)*imscale;

% number of tiles in X and Y
nX = round(size(image, 2)/unitsize);
nY = round(size(image, 1)/unitsize);

% center positions
gridX = round((meshgrid(1:nX, 1:nY)-1)*unitsize + tilesize*imscale*.5);
gridY = round((meshgrid(1:nY, 1:nX)-1)*unitsize + tilesize*imscale*.5)';
gridpos = [gridX(:), gridY(:)]/imscale;

[~, ~, imOrder] = imfiles_in_grid('', '', '', 0, [nX nY], imorder);

% empty tiles
for i = 1:nX*nY    
    if image(gridY(i),gridX(i),1) < 2
        imOrder(i) = 0;
    end
end

[uTiles, iOrder] = unique(imOrder);
if uTiles(1) == 0
    iOrder = iOrder(2:end);
end

% tile position
tileX = meshgrid(1:nX, 1:nY);
tileY = meshgrid(1:nY, 1:nX)';
tilePosYX = [tileY(iOrder), tileX(iOrder)];

% grid position
gridpos = gridpos(iOrder,:);

end

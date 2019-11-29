function [Coord, Coord_write, blank] = drawpolygon(I, nROIs, scale, isYInverted)
% [Coord, Coord_write, blank] = drawpolygon(I, nROIs, scale, isYInverted)
% draw polygons on an image and get coordinates
% Xiaoyan, 2018

Coord = {}; Coord_write = []; blank = [];

if nargin <= 2 
    scale = 1;
end
    
if nargin < 4
    isYInverted = 0;
end

for i = 1:nROIs
    figure;
    set(gcf, 'units', 'normalized', 'position', [0.01 0.05 0.98 0.85]);    
    imshow(I, []);
    if isYInverted
        set(gca, 'YDir', 'normal')
    end
    
    % show previously drawn polygons
    hold on;
    for j = 1:i-1
        plot(Coord{2,j}*scale, Coord{3,j}*scale, 'linewidth', 2);
    end
    
%         h = imfreehand;
    h = impoly;

    % wait until double click to confirm polygon
    position = wait(h);
    
    position = floor(position/scale);
    Coord = [Coord,...
        [{['ROI' num2str(i)]};...
        {[position(:,1);position(1,1)]'};...
        {[position(:,2);position(1,2)]'}]];
    
    Coord_write = [Coord_write; zeros(size(position,1),1)+i,position];
    
    % if no more than three points found in polygon, add to blank list
    if size(position,1) <= 2
        blank = [blank; i];
    end
    close(gcf);
    
end

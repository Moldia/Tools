function [Coord, Coord_write, blank] = drawpolygon_density(density, nROIs, imsize)
% [Coord, Coord_write, blank] = drawpolygon_density(density, nROIs, imsize)
% draw polygons from an estimated density plot and get coordinates
% coordinates match the scaled background image
% Xiaoyan, 2018

Coord = {}; Coord_write = []; blank = [];


for i = 1:nROIs
    figure;
    set(gcf, 'Units', 'Normalized', 'Position', [0.01 0.05 0.98 0.85]);
    imshow(density, []);
    colormap(gca, parula);

    %     h = imfreehand;
    h = impoly(gca);
    pos = getPosition(h);
    % wait until double click to confirm polygon
    position = wait(h);
    position = position*5;
    
    Coord = [Coord,...
        [{['ROI' num2str(i)]};...
        {[position(:,1);position(1,1)]'};...
        {[position(:,2);position(1,2)]'}]];
    
    Coord_write = [Coord_write; ...
        zeros(size(position,1),1)+i,position];
    
    if size(position,1)==1 || size(position,1)==2
        blank = [blank;i];
    end
    
    close(gcf);
end

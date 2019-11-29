% visualize pooled tSNE results on top of background images
% Xiaoyan, 2018

clear; close all;

%% modify here
tSNE_results = 'tSNE_3D.csv';
hexbin_counts = 'E:\Whole_organoid_pseudoAnchor_v5\Decoding\pooled_bincounts_count.csv';   % already binned data, with the row names specifying sample and bin#
hexbin_position = 'E:\Whole_organoid_pseudoAnchor_v5\Decoding\pooled_bincounts_binpos.csv';
hexbin_size = 200;	% if single-cell data, hexbin_size = 0;
images_to_use = 'background_images.txt';     % txt file, one image file per row, follow exactly the original order in files_to_pool file
scale = 1;  % image scale, applies to all images

%% do not modify

% get sample and bin names
tableCount = readtable(hexbin_counts, 'ReadVariableNames', 1);
[uSamples, ~, iSample] = unique(cellfun(@(v) v(1:strfind(v, '_hexbin')-1), table2cell(tableCount(:,1)), 'uni', 0));

% get position
pos = importdata(hexbin_position);
pos = pos.data;

% image files
fid = fopen(images_to_use, 'r');
imfiles = textscan(fid, '%s', 'delimiter', '\n');
imfiles = imfiles{:};
fclose(fid);

% load tSNE results
Y = csvread(tSNE_results);

% visualize tSNE in RGB
if ~hexbin_size;	hexbin_size = 5;    end
Yrgb = rgbscale(Y);

for s = 1:numel(uSamples)
    figure; 
    try
        imshow(imfiles{s});
    end
    hold on;
    for i = 1:nnz(iSample==s)
        pos_sample = pos(iSample==s,:);
        Yrgb_sample = Yrgb(iSample==s,:);
        [gridR, gridL, xypos] = hexbin_coord2grid(pos_sample(i,:), hexbin_size);
        [vy, vx] = dot2poly(pos_sample(i,2)*scale,pos_sample(i,1)*scale,...
            hexbin_size*scale, 6);
        patch(vx, vy, Yrgb_sample(i,:), 'facealpha', .2);
    end
    title(uSamples{s});
    axis image
    axis off
    drawnow;
end

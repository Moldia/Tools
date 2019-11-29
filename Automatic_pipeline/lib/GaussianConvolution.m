%% density plot (guassian smoothing)
%  Xiaoyan, 2015-4-7

clear;
% close all;
drawnow;

%% parameters
decoded_file = 'K:\161230_161220KI_3-1\Kenneth\hippocampi\hippocampus1\spots_ROI1.csv';
image = 'K:\161230_161220KI_3-1\Kenneth\hippocampi\hippocampus1\Ab_c1_ROI1.tif';    % important for size
name_density = 'Pvalb';
bandwid = 50;

%% transcripts
[name, pos] = getinsitudata(decoded_file);

% unique transcripts
[uniName, ~, idxName] = unique(name);
[p, q] = hist(idxName,unique(idxName));

%% image size
img = imfinfo(image);
imsize = [img.Height, img.Width];

%% gaussian smoothing
idx = find(strcmp(uniName, name_density));
if isempty(idx)
    error('No specified transcript detected in the input file');
end
pos_gaussian = pos(idxName==idx, :);

temp = floor(pos_gaussian/5);
temp(temp==0) = 1;
Itemp = accumarray(fliplr(temp), 1, floor(imsize/5));
fh = fspecial('gaussian', bandwid*2, bandwid/5);
Itemp = imfilter(Itemp, fh);

figure;
imshow(Itemp/max(fh(:)),[]);
colormap(hot);
colorbar
title(name_density);
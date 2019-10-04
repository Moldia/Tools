%% density plot (guassian smoothing)
%  Xiaoyan, 2015-4-7

clear;
% close all;
drawnow;

%% parameters
decoded_file = 'E:\PROOOJECTS\test_dataset\QT_0.35_0.004_details.csv';
image = 'C:\Worky\Temp\TissueA_Fluo_staining.png';   % important for size
name_density = 'HER2';
bandwid = 50;

%% transcripts
[name, pos] = getinsitudata(decoded_file);

% unique transcripts
[uNames, ~, iName] = unique(name);
% [p, q] = hist(iName, unique(iName));

%% image size
img = imfinfo(image);
imsize = [img.Height, img.Width];

%% gaussian smoothing
idx = find(strcmp(uNames, name_density));
if isempty(idx)
    error('No specified transcript detected in the input file');
end
pos_gaussian = pos(iName==idx, :);

temp = floor(pos_gaussian/5);
temp(temp==0) = 1;
Itemp = accumarray(fliplr(temp), 1, floor(imsize/5));
Imask = imdilate(Itemp, strel('disk', 4));
fh = fspecial('gaussian', bandwid*2, bandwid/5);
Igauss = imfilter(Itemp, fh);

figure;
imshow(Igauss/max(fh(:)).*Imask, []);
cmap = colormap(jet);
cmap(1,:) = 1;
colormap(gca, cmap);
colorbar
title(name_density);
drawnow;

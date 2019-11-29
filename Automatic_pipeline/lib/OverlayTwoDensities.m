% overlay density estimates of two genes
% Xiaoyan, 2017

clear;
close all;
drawnow;

%% modify here
decoded_file = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\1\QT_0.6_details_noNNNN.csv';
image = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\1\Sample 1_Base 3_resize20_c1.jpg';
scale = .2;
genes_density = {'Homo sapiens microtubule associated protein 2 (MAP2)', 'Homo sapiens cyclin dependent kinase inhibitor 1A (CDKN1A)','Homo sapiens nucleolar and spindle associated protein 1 (NUSAP1)'};       % first in red, second in cyan
bandwidth = 100;  % in original scale
%% do not modify

% load
[name, pos] = getinsitudata(decoded_file);

% unique transcripts
[uNames, ~, iName] = unique(name);
[p, q] = hist(iName, 1:length(uNames));

% image size
img = imfinfo(image);
imsize = [img.Height, img.Width];

layer = [1 0 0; 0 1 1];

% densities
I = zeros(ceil(imsize(1)/5), ceil(imsize(2)/5), 3, 'double');
for i = 1:2
    density = gene_kde(name, pos, genes_density{i}, bandwidth, image, scale);
    density = (density - min(density(:)))/max(density(:));   
    I = I + cat(3, density*layer(i,1), density*layer(i,2), density*layer(i,3));
end

figure;
imshow(I/max(I(:))*2);
title([genes_density{1} ' - ' genes_density{2} ' densities' ]);

% absolute counts
I = zeros(floor(imsize(1)/5), floor(imsize(2)/5), 3, 'double');
for i = 1:2
    pos_density = pos(iName==find(strcmp(genes_density{i}, uNames)),:);
    temp = floor(pos_density*scale/5);
    temp(temp==0) = 1;
    Itemp = accumarray(fliplr(temp), 1, floor(imsize/5));
    fh = fspecial('gaussian', bandwidth*2, bandwidth/5);
    smoothed = imfilter(Itemp, fh);
    smoothed = smoothed/max(fh(:));     % normalize to peak height of gaussian filter
    I = I + cat(3, smoothed*layer(i,1), smoothed*layer(i,2), smoothed*layer(i,3));
end

figure;
imshow(I/max(I(:)));
title([genes_density{1} ' - ' genes_density{2} ' gaussian smoothed' ]);

% kernel density estimation of a specific gene
% Xiaoyan, 2018

clear;
close all;
drawnow;

%% parameters
decoded_file = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\1\QT_0.6_details_noNNNN.csv';
image = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\1\Sample 1_Base 3_resize20_c1.jpg';   % important for size
scale = .2;      % image scale
gene_density = 'Homo sapiens tubulin beta 3 class III (TUBB3)';
bandwid = 100;   % in original scale

%%
% all transcripts
[name, pos] = getinsitudata(decoded_file);

% density estimation plot
density = gene_kde(name, pos, gene_density, bandwid, image, scale);

figure;
imshow(density, []);
colormap(gca, parula);
title(gene_density);

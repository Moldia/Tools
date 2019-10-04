% kernel density estimation of a specific gene
% Xiaoyan, 2018

clear;
close all;
drawnow;

%% parameters
decoded_file = 'E:\PROOOJECTS\test_dataset\QT_0.35_0.004_details.csv';
image = 'E:\PROOOJECTS\test_dataset\860502_1_align.png';   % important for size
scale = .2;      % image scale
gene_density = 'ACTB';
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

% plotting on top of the density estimate
% Xiaoyan, 2018

close all; drawnow;

%% modify here
decoded_file = 'E:\PROOOJECTS\test_dataset\QT_0.35_0.004_details.csv';
image = 'E:\PROOOJECTS\test_dataset\860502_1_align.png';    % important for size
scale = .2;      % image scale

gene_density = 'COL3A1';
genes_to_plot = {'VIM', 'HER2', 'TP53', 'CD68', 'MUC1'};      % if empty, all reads will be plotted

bandwid = 200;   % in original scale

%% do not modify

% load
[name, pos] = getinsitudata(decoded_file);

% kde
density = gene_kde(name, pos, gene_density, bandwid, image, scale);

% plot
figure;
imshow(density, []);
colormap(gca, parula);
hold on;
plotall(name, correctcoord(pos, scale/5));
update_legend(gca, genes_to_plot);
title([gene_density ' density']);

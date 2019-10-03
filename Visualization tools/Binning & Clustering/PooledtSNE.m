% use bin count data to run PCA and tSNE
% Xiaoyan, 2018

clear; close all;


%% modify here
hexbin_counts = 'pooled_bincounts_count.csv';  % already binned data (pooled or not), requires header, compatible with output from BatchHexBin
hexbin_position = 'pooled_bincounts_binpos.csv';    % bin position file generated at the same time as bin count file
hexbin_size = 500;   % in piexels, if single-cell data use hexbin_size = 0
output_directory = '';

%% do not modify

% import data
tableCount = readtable(hexbin_counts, 'ReadVariableNames', 1);

% original column names
cNames = tableCount.Properties.VariableNames;

% make sure there are no multiple entries of the same gene
assert(numel(cNames)== numel(unique(cNames)),...
    'Column names are not unique!')

% create checkboxes and get selected values
cbValues = checkboxes(cNames);
idx = cellfun(@(v) find(strcmp(v, cNames)), cbValues(:,1));
isSelected = false(numel(idx), 1);
isSelected(idx) = cell2mat(cbValues(:,2));
cGenes = table2array(tableCount(:,isSelected));
genes = cNames(isSelected)'

% PCA
[coeff, score, latent] = pca(cGenes);
% visualize first two components
figure, biplot(coeff(:,1:2), 'Scores', score(:,1:2), 'VarLabels', genes);
title('top two principle components');

% tSNE in MATLAB
% ONLY >=R2018a
seeds = 1e-4*randn(size(cGenes,1), 3);
Y = tsne(cGenes, 'NumDimensions', 3, 'NumPCAComponents', 50, 'Perplexity', 30,...
    'Standardize', 1, 'LearnRate', 1000, 'Verbose', 1, 'InitialY', seeds); 
[uSamples, ~, iSample] = unique(cellfun(@(v) v(1:strfind(v, '_hexbin')-1), table2cell(tableCount(:,1)), 'uni', 0));
figure, gscatter(Y(:,1),Y(:,2), iSample);
legend(uSamples);
title({'tSNE dim reduction to three, shown first two', 'color-coded by samples'});

% get position
pos = importdata(hexbin_position);
pos = pos.data;

% visualize tSNE in RGB (no background)
if ~hexbin_size;	hexbin_size = 5;    end
Yrgb = rgbscale(Y);
for s = 1:numel(uSamples)
    figure;
    hold on;
    scale = 1;
    for i = 1:nnz(iSample==s)
        pos_sample = pos(iSample==s,:);
        Yrgb_sample = Yrgb(iSample==s,:);
        [gridR, gridL, xypos] = hexbin_coord2grid(pos_sample(i,:), hexbin_size);
        [vy, vx] = dot2poly(pos_sample(i,2)*scale,pos_sample(i,1)*scale,...
            hexbin_size*scale, 6);
        patch(vx, vy, Yrgb_sample(i,:));
        %     plot(pos(idx,1), pos(idx,2), 'o',...
        %         'color', Yrgb(i,:),...
        %         'linewidth', 3);
    end
    title(uSamples{s});
    axis image
    axis off
    drawnow;
end

% write
mkdir(output_directory);
csvwrite(fullfile(output_directory, 'tSNE_3D.csv'), Y);
csvwrite(fullfile(output_directory, 'tSNE_initial.csv'), seeds);


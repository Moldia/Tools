% unsupervised clustering (kmeans) of spatial data
% transcript counts in every hexbin is normalized by the maximum counts
% Xiaoyan, 2017

%% modify here
hexagon_radius = 150;   % in pixel
decoded_file = 'E:\Organoid 1\Organoid 1S\QT_0.45_1e-07_details.csv';
num_clusters = 3;
output_directory = '';

background_image = '';    % can be empty if not needed for visualization
scale = 0.2;    % image scale

%% do not modify

% transcripts
[name, pos] = getinsitudata(decoded_file);
[name, pos] = removereads(name, 'NNNN', pos);
[uNames, ~, iName] = unique(name);

% reads in hexagon
[inbin, bincenters, vx, vy] = hexbin(pos, hexagon_radius);
nonEmptyGrids = unique(inbin);

% gene counts
cGenes_nonzero = zeros(length(uNames), size(bincenters,1));
for i = 1:numel(nonEmptyGrids)
    cGenes_nonzero(:,i) = (hist(iName(inbin==nonEmptyGrids(i)), 1:numel(uNames)))';
end

% write files (all data)
cNames = {'grid_num', uNames{:}, 'center_x', 'center_y'};

fid = fopen(fullfile(output_directory, 'HexbinClustering_GeneCount_all.csv'), 'w');
fprintf(fid, lineformat('%s', numel(cNames)), cNames{:});
towrite = num2cell([nonEmptyGrids, cGenes_nonzero', bincenters]');
fprintf(fid, lineformat('%d', numel(cNames)), towrite{:});
fclose(fid);

% select genes for clustering
cbValues = checkboxes(uNames);
idx = cellfun(@(v) find(strcmp(v, uNames)), cbValues(:,1));
isSelected = false(numel(idx), 1);
isSelected(idx) = cell2mat(cbValues(:,2));

cGenes_selected = cGenes_nonzero(isSelected,:);
emptygrid = sum(cGenes_selected,1)==0;
nonEmptyGrids = nonEmptyGrids(~emptygrid);
cGenes_selected = cGenes_selected(:,~emptygrid);

% normalized by max
cGenes_maxnorm = bsxfun(@rdivide, cGenes_selected, max(cGenes_selected,[],2));

% normalization by sum
cGenes_sumnorm = bsxfun(@rdivide, cGenes_selected, sum(cGenes_selected,1));


% kmeans
disp('Starting kmeans clustering with 500 replicates..');
[iCluster, centroid] = kmeans(cGenes_maxnorm', num_clusters,...
    'Distance', 'sqeuclidean', 'Replicates', 500);

% write files (only selected)
cNames = {'grid_num', uNames{isSelected}, 'center_x', 'center_y'};

fid = fopen(fullfile(output_directory, 'HexbinClustering_GeneCount.csv'), 'w');
fprintf(fid, lineformat('%s', numel(cNames)), cNames{:});
towrite = num2cell([nonEmptyGrids, cGenes_selected', bincenters(~emptygrid,:)]');
fprintf(fid, lineformat('%d', numel(cNames)), towrite{:});
fclose(fid);

fid = fopen(fullfile(output_directory, 'HexbinClustering_GeneCount_SumNorm.csv'), 'w');
fprintf(fid, lineformat('%s', numel(cNames)), cNames{:});
towrite = num2cell([nonEmptyGrids, cGenes_sumnorm', bincenters(~emptygrid,:)]');
fprintf(fid, lineformat('%d', numel(cNames)), towrite{:});
fclose(fid);

cNames = {'grid_num', uNames{isSelected}, 'center_x', 'center_y', 'cluseter_id'};
fid = fopen(fullfile(output_directory, 'HexbinClustering_GeneCount_MaxNorm.csv'), 'w');
fprintf(fid, lineformat('%s', numel(cNames)), cNames{:});
towrite = num2cell([nonEmptyGrids, cGenes_maxnorm', bincenters(~emptygrid,:), iCluster]');
fprintf(fid, lineformat('%d', numel(cNames)), towrite{:});
fclose(fid);

% visualization
vx = vx(:,~emptygrid);
vy = vy(:,~emptygrid);
figure; 
try
    image = imread(background_image);
    imshow(image);
catch
    axis image
    set(gca, 'YDir', 'reverse');
end
col = {'red' 'green' 'blue' 'yellow' 'cyan'};
hold on;
for k = 1:max(iCluster)
    bins = find(iCluster==k);
    for i = 1:length(bins)
        fill(vx(:,bins(i))*scale, vy(:,bins(i))*scale, col{k}, 'facealpha', 0.3);
    end
end
drawnow;

% heatmap of normalized counts and barplot of cluster centroid position
figure;
ax1 = subplot(121);
[~, idxSort] = sort(iCluster);
[~, idxFirst] = unique(iCluster(idxSort));
imagesc(cGenes_maxnorm(:,idxSort));
hold on;
plot(repmat(idxFirst', 2, 1), repmat([0; nnz(isSelected)+1], 1, numel(idxFirst)),...
    'r', 'linewidth', 2);
xlabel('bin');
set(gca, 'ytick', 1:nnz(isSelected), 'yticklabel', uNames(isSelected), 'fontsize', 5);
title('normailzed bin count data')
colorbar

ax2 = subplot(122);
bh = barh(centroid');
set(gca, 'ytick', 1:nnz(isSelected), 'yticklabel', uNames(isSelected),...
    'ylim',[0 nnz(isSelected)+1], 'fontsize', 5, 'ydir', 'reverse');
for i = 1:length(bh)
    set(bh(i), 'facecolor', rgb(col{i}), 'edgecolor', [.2 .2 .2], 'linewidth', .1);
end
box off
xlabel('normalized gene count')
title('centroid location of clusters')
legend(catstrnum('Cluster ', 1:max(iCluster)))

linkaxes([ax1, ax2], 'y');

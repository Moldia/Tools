% cluster already processed hexbin count file
% Xiaoyan, 2018


%% modify here
hexbin_count_file = 'J:\NextPipeline\HexbinClustering_GeneCount_all.csv';
num_clusters = 3;
output_directory = '';
background_image = '';
scale = 1;      % scale of background image

%% do not modify

% import
data = importdata(hexbin_count_file);
if ismember('cluseter_id', data.colheaders)
    cGenes_selected = data.data(:,2:end-3)';
    binpos = data.data(:,end-2:end-1);
    
    cGenes_maxnorm = cGenes_selected;
    cNames = data.colheaders;
    
else
    cGenes_selected = data.data(:,2:end-2)';
    binpos = data.data(:,end-1:end);
    
    % normalized by max
    cGenes_maxnorm = bsxfun(@rdivide, cGenes_selected, max(cGenes_selected,[],2));
    
    cNames = [data.colheaders, {'cluseter_id'}];
    
end

% estimate radius
radius = unique(binpos(:,2));
radius = min(abs(diff(radius)))*2/3;

[vy, vx] = dot2poly(binpos(:,2), binpos(:,1), radius, 6);


% kmeans
disp('Starting kmeans clustering with 100 replicates..');
[iCluster, centroid] = kmeans(cGenes_maxnorm', num_clusters,...
    'Distance', 'sqeuclidean', 'Replicates', 100);

% write files
fid = fopen(fullfile(output_directory, 'HexbinClustering_GeneCount_MaxNorm_SecondaryClustering.csv'), 'w');
fprintf(fid, lineformat('%s', numel(cNames)), cNames{:});
towrite = num2cell([data.data(:,1), cGenes_maxnorm', binpos, iCluster]');
fprintf(fid, lineformat('%d', numel(cNames)), towrite{:});
fclose(fid);

% visualization
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
plot(repmat(idxFirst', 2, 1), repmat([0; size(cGenes_maxnorm,1)+1], 1, numel(idxFirst)),...
    'r', 'linewidth', 2);
xlabel('bin');
set(gca, 'ytick', 1:size(cGenes_maxnorm,1), 'yticklabel', cNames(2:end-3), 'fontsize', 5);
title('normailzed bin count data')
colorbar

ax2 = subplot(122);
bh = barh(centroid');
set(gca, 'ytick', 1:size(cGenes_maxnorm,1), 'yticklabel', cNames(2:end-3),...
    'ylim',[0 size(cGenes_maxnorm,1)+1], 'fontsize', 5, 'ydir', 'reverse');
for i = 1:length(bh)
    set(bh(i), 'facecolor', rgb(col{i}), 'edgecolor', [.2 .2 .2], 'linewidth', .1);
end
box off
xlabel('normalized gene count')
title('centroid location of clusters')
legend(catstrnum('Cluster ', 1:max(iCluster)))

linkaxes([ax1, ax2], 'y');

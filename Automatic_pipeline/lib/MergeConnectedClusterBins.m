% merge spatially connected hexagonal bins (requires clustered data)
% Xiaoyan, 2017

clear; close all;

%% modify here
hexbin_cluster_file = 'HexbinClustering_GeneCount_MaxNorm.csv';    % last three columns: x pos, y pos, cluster id
hexagon_radius = 600;   % in pixels
output_directory = '';
min_cluster_member = 20;     % minimum number of bins in a connected spatial patch

background_image = '';    % can be empty if not needed for visualization
scale = 0.2;    % image scale

%% do not modify

% load data
data = importdata(hexbin_cluster_file);
iCluster = data.data(:,end);
hexcentroid = data.data(:,end-2:end-1);
[gridR, gridL, xypos] = hexbin_coord2grid(hexcentroid, hexagon_radius);

% connected bins
connected = zeros(length(xypos), max(iCluster)+1, 'uint8');
for c = 1:max(iCluster)
    clear labelPixelR labelPixelL

    % connect left, right, upper right, lower right
    xyposR = ceil(xypos);
    [~, sortR] = sortrows(xyposR);    
    gridSub = zeros(size(gridR), 'uint8');
    nonEmpty = find(gridR(:));
    gridSub(nonEmpty(iCluster(sortR)==c)) = 1;
    cc = bwconncomp(gridSub, 4);
    labelR = labelmatrix(cc);
    labelPixelR(sortR) = labelR(nonEmpty);
     
    % connect left, right, upper left, lower left
    xyposL = floor(xypos);
    [~, sortL] = sortrows(xyposL);
    gridSub = zeros(size(gridL), 'uint8');
    nonEmpty = find(gridL(:));
    gridSub(nonEmpty(iCluster(sortL)==c)) = 1;
    cc = bwconncomp(gridSub, 4);
    labelL = labelmatrix(cc);
    labelPixelL(sortL) = labelL(nonEmpty);

    % merge two connected maps
    toMerge = unique([labelPixelR(:), labelPixelL(:)], 'rows');
    toMerge = toMerge(all(toMerge,2),:);
    merged = groupPQ([(1:size(toMerge,1))', toMerge]);
    
    labelPixel = zeros(length(xypos), 1, 'uint8');
    for i = 1:length(merged)
        labelPixel(labelPixelR==toMerge(i,1)) = merged(i);
        labelPixel(labelPixelL==toMerge(i,2)) = merged(i);
    end
    
    % remove patchs with <min_cluster_member bins
    [uLabels, ~, iLabel] = unique(labelPixel);
    cLabel = hist(iLabel, 1:numel(uLabels));
    labelPixel(ismember(iLabel, find(cLabel<min_cluster_member))) = 0;
    [~, ~, labelPixel] = unique(labelPixel);
    if ~uLabels(1)
        labelPixel = labelPixel - 1;
    end
    
    connected(labelPixel~=0,1) = c;
    connected(:,c+1) = labelPixel;
end

% visualize all clusters
figure; 
col = {'red' 'green' 'blue' 'yellow' 'cyan'};
subplot(231);
try
    image = imread(background_image);
    imshow(image);
catch
    axis image
    set(gca, 'YDir', 'reverse');
end
hold on;
[vy, vx] = dot2poly(hexcentroid(:,2)*scale,hexcentroid(:,1)*scale,...
    hexagon_radius*scale, 6);

for c = 1:max(iCluster)
    fill(vx(:,iCluster==c), vy(:,iCluster==c), col{c},...
        'facealpha', 0.3, 'edgecolor', ' none');
end
xvrange = get(gca, 'xlim');
yvrange = get(gca, 'ylim');
title('all clusters');

% get polygon outlines after merging and plot
Outlines = cell(max(iCluster),1);
linecol = lines(max(reshape(connected(:,2:end), [], 1)));

for c = 1:max(iCluster)
    subplot(2,3,c+1);
    try
        imshow(image);
    catch
        axis image
        set(gca, 'YDir', 'reverse',...
            'xlim', xvrange,...
            'ylim', yvrange);
    end
    hold on;
    title(['Cluster' num2str(c) ' spatial subclusters']);
    
    for l = unique(connected(:,c+1))'
        if l
            hexpossub = hexcentroid(connected(:,c+1)==l,:);
            
            hold on;
            [vy, vx] = dot2poly(hexpossub(:,2)*scale, hexpossub(:,1)*scale,...
                hexagon_radius*scale, 6);
            fill(vx, vy, col{c}, 'facealpha', .1, 'edgecolor', 'none');
            
            outline = [];
            for p = 1:size(hexpossub, 1)
                [vy, vx] = dot2poly(hexpossub(p,2), hexpossub(p,1), hexagon_radius, 6);
                pairx = [vx(1:end-1), vx(2:end)];
                pairy = [vy(1:end-1), vy(2:end)];
                
                outline = [outline; pairx, pairy];
            end
            
            % somehow the original numbers do not have the same precision
            % and this causes problem in comparing numbers
            contour_rounded = round(outline, -floor(log10(hexagon_radius)));
          
            for i = 1:size(contour_rounded,1)
                pairx = contour_rounded(i,1:2);
                pairy = contour_rounded(i,3:4);
                
                [~, idxY] = sort(pairy, 2);
                [~, idxX] = sort(pairx, 2);
                
                if isequal(idxX, idxY)
                    outline(i,:) = outline(i,[idxX, idxY+2]);
                    contour_rounded(i,:) = [pairx(idxX), pairy(idxY)];
                elseif pairx(1)==pairx(2)
                    outline(i,:) = outline(i,[idxY, idxY+2]);
                    contour_rounded(i,:) = [pairx(idxY), pairy(idxY)];
                else
                    outline(i,:) = outline(i,[idxX, idxX+2]);
                    contour_rounded(i,:) = [pairx(idxX), pairy(idxX)];
                end
            end
            
            [uPairs, iFirst, iPair] = unique(contour_rounded, 'rows');
            cPairs = hist(iPair, 1:size(uPairs,1));
            outline = outline(iFirst(cPairs==1,:),:);

            plot(outline(:,1:2)'*scale, outline(:,3:4)'*scale,...
                'color', linecol(l,:), 'linewidth', 1);
            
            Outlines{c} = [Outlines{c}, {outline}];
            
        end
        drawnow;
    end
end

save(fullfile(output_directory, 'MergedContours.mat'), 'Outlines');

% specify a gene to create density estimate
% threshold the desity estimate to get binary mask
% extract reads within the mask
% Xiaoyan, 2018

clear;
close all;
drawnow;

%% parameters
decoded_file = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\3\QT_0.65_details_noNNNN.csv';
background_image = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\3\Sample 3_Base 3_resize20_c1.jpg';     % important for size
scale = .2;  % image scale

name_density = 'Homo sapiens kinesin family member 23 (KIF23)';   
name_exclude = ''; % exclude some genes from the plot and output files
bandwidth = 300;
density_threshold = .30;   % between 0 and 1, if 1, automatic thresholding
area_threshold = 200;     % minimal area for a segmented region, in original scale

show_image = 1; % show background image for plotting AND show ROI on top of the image
output_directory = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\1';

%%
% load
[name, pos, tilepos, general, qscore] = ...
    getinsitudata(decoded_file, 2, 1, 3, 4, 5);

% remove specified transcripts
[name, pos, tilepos, general, qscore] = ...
    removereads(name, name_exclude, pos, tilepos, general, qscore);
[uNames, ~, iName] = unique(name);
cName = hist(iName, 1:length(uNames));
[~, idxSort] = sort(cName, 'descend');

posScaled = correctcoord(pos, scale);

% kde
density = gene_kde(name, pos, name_density, bandwidth, background_image, scale);
densityScaled = density/max(density(:));

% binarize density
if density_threshold == 1
    densityBW = im2bw(densityScaled, graythresh(densityScaled));
    fprintf('Threshold used: %.3f\n', graythresh(densityScaled));
else
    densityBW = im2bw(densityScaled, density_threshold);
end

% segment objects and threshold based on size
[densityLabel, clusterCentroid, clusterProps] =...
    filterobjects(densityBW, 'Area', area_threshold*.2);
clusterBoundary = bwboundaries(densityLabel);

% original scale
clusterCentroid = clusterCentroid*5/scale;
clusterProps = clusterProps*(5/scale)*(5/scale);
clusterBoundary = cellfun(@(v) v*5/scale, clusterBoundary, 'uni', 0);

% ROI visualization
col = lines(length(clusterBoundary));
if show_image
    disp('loading image...')
    im = imread(background_image);
    figure; imshow(im, []); hold on;
else
    plotonblank;
end
for i = 1:length(clusterBoundary)
    plot(clusterBoundary{i}(:,2)*scale, clusterBoundary{i}(:,1)*scale,...
        'color', col(i,:), 'linewidth', 2);
end
text(clusterCentroid(:,1)*scale, clusterCentroid(:,2)*scale,...
    cellfun(@num2str, num2cell(1:max(densityLabel(:))), 'uni', 0)',...
    'fontsize', 12, 'color', 'red', 'horizontalalign', 'center');

% find reads in ROI
inROI = readsinroi(posScaled*.2, densityLabel);

% plot reads
figure;
plotall(name(logical(inROI)), pos(logical(inROI),:), background_image, scale);
drawnow;

% ROI (hole-filled) coordinates
filledClusterBoundary = bwboundaries(imfill(densityLabel));
filledClusterBoundary = cellfun(@(v) v*5/scale, filledClusterBoundary,...
    'uni', 0);
boundarylength = [0; cellfun(@length, filledClusterBoundary)];
towrite = zeros(sum(boundarylength),3);
izero = 0;
for i = 1:length(filledClusterBoundary)
    izero = izero+boundarylength(i);
    towrite(izero+1:izero+boundarylength(i+1),1) = i;
    towrite(izero+1:izero+boundarylength(i+1),2:3) = fliplr(filledClusterBoundary{i});
end

% write ROI coordinates
mkdir(output_directory);
fid = fopen(fullfile(output_directory, 'DensityThreshCoordinates.csv'),'w');
fprintf(fid,'Polygon id,x coordiates,y coordinates\n');
fprintf(fid,'%d,%d,%d\n',towrite');
fclose(fid);

% counts
histReads = hist3(...
    [iName, inROI],...
    [{1:length(uNames)}, {0:max(densityLabel(:))}]);
histReads = histReads(idxSort,:);
towrite = [uNames(idxSort), num2cell(histReads(:,2:end))]';
fid = fopen(fullfile(output_directory, 'DensityThreshCounts.csv'), 'w');
fmt = lineformat('%s', max(densityLabel(:)));
headerline = catstrnum('ROI', 1:max(densityLabel(:)));
fprintf(fid, ['gene,', fmt], headerline{:});
fmt = lineformat('%d', max(densityLabel(:)));
fprintf(fid, ['%s,', fmt], towrite{:});
fprintf(fid,'\nEstimated ROI area,');
fprintf(fid, fmt, clusterProps);
fclose(fid);

% details file with ROI number attached
towrite = [name, name,...     % to keep same format as input (for convenience...)
    num2cell([pos, tilepos, general, qscore, inROI])];
% towrite = towrite(inROI~=0,:)';
towrite = towrite';
fid = fopen(fullfile(output_directory, 'DensityThreshDetails.csv'), 'w');
fprintf(fid, 'Name,Name,GlobalX,GlobalY,TilePos,GeneralStain,Qscore,inROI\n');
fprintf(fid, ['%s,%s,', lineformat('%d', 6)], towrite{:});
fclose(fid);

% barplot
figure; 
bh = bar3(histReads(:,2:end)./repmat(clusterProps',length(uNames),1), 1);
for i = 1:length(bh)
    bh(i).CData = bh(i).ZData;
end
view(2)
axis normal
set(gca, 'YTick', 1:length(uNames), 'YTickLabel', uNames(idxSort),...
    'YLim', [.5 length(uNames)+.5], 'XTick', 1:2:max(densityLabel(:)),...
    'FontSize', 5);
xlabel('ROI #');
title('Counts within ROI/ROI area');
set(gcf, 'Position', [200 200 800 600]);

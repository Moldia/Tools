% draw polygons to create ROI (on the estimated density plot)
% polygons can be adjusted after they are initially created
% double click to confirm a polygon
% Xiaoyan, 2018

clear;
close all;

%% modify here
decoded_file = 'E:\PROOOJECTS\test_dataset\QT_0.35_0.004_details.csv';
background_image = 'E:\PROOOJECTS\test_dataset\860502_1_align.png';     % important for size
scale = .2;  % image scale

roi_number = 3;   % how many ROI regions you want to extract

gene_density = 'COL3A1';   
genes_to_exclude = {'NNNN'};    % exclude from output files
bandwidth = 200;

show_image = 1; % 0-show ROIs on density; 1-show ROIs on original image
output_directory = 'E:\PROOOJECTS\test_dataset\ROI_DrawDensity';

%% do not modify

% specify colors
col = {'white' 'yellow' 'blue' 'tomato' 'palevioletred' 'bisque' 'pink' 'mediumspringgreen'};

% load
[name, pos] = getinsitudata(decoded_file);

% remove specified transcripts
[name, pos] = removereads(name, genes_to_exclude, pos);
[uNames, ~, iName] = unique(name);
cName = hist(iName, 1:length(uNames));
[~, idxSort] = sort(cName, 'descend');

posScaled = correctcoord(pos, scale);

% kde
density = gene_kde(name, pos, gene_density, bandwidth, background_image, scale);

% create ROIs
[Coord, Coord_write, blank] = drawpolygon_density(density, roi_number);
Coord_write(:,2:3) = Coord_write(:,2:3)/scale;

% remove empty polygons
if ~isempty(blank)
    Coord(:,blank) = [];
end
ROI_exist = 1:roi_number;
ROI_exist(blank) = [];

% ROI visualization
if length(col)<roi_number
    col = repmat(col, 1, ceil(roi_number/length(col)));
end
if show_image
    disp('loading image...')
    im = imread(background_image);
    figure; imshow(im, []); hold on;
    for i = 1:length(ROI_exist)
        plot(Coord{2,i}, Coord{3,i}, 'LineWidth', 2, 'Color', rgb(col{ROI_exist(i)}));
    end
else
    figure; imshow(density, []); 
    colormap parula; hold on;
    for i = 1:length(ROI_exist)
        plot(Coord{2,i}/5, Coord{3,i}/5, 'LineWidth', 2, 'Color', rgb(col{ROI_exist(i)}));
    end
    title(gene_density);
end
legend(Coord(1,:), 'Location', 'NorthEastOutside', 'Color',[.6 .6 .6]);
drawnow;

% write coordinates to file
mkdir(output_directory);
fid = fopen(fullfile(output_directory, 'DensityROICoordinates.csv'),'w');
fprintf(fid,'Polygon id,x coordiates,y coordinates\n');
fprintf(fid,'%d,%d,%d\n',Coord_write');
fclose(fid);

% count transcripts within polygons
[ROI_count, ROI_freq, ROI_proportion, ROI_area, In] = ...
    count_in_polygon(length(ROI_exist), Coord,...
    posScaled, uNames, iName, cName);

% heatmap
% figure;
% bh = bar3(ROI_freq(idxSort,:), 1);
% for i = 1:length(bh)
%     bh(i).CData = bh(i).ZData;
% end
% view(2)
% set(gca, 'YTick', 1:length(uNames), 'YTickLabel', uNames(idxSort),...
%     'YLim', [.5 length(uNames)+.5], 'XTick', ROI_exist,...
%     'FontSize', 5);
% axis normal
% xlabel('ROI #');
% title('Normalized counts');
% set(gcf, 'Position', [200 200 800 600]);

% bar plot
figure;
bh = bar(ROI_freq);
for i = 1:length(bh)
    set(bh(i), 'FaceColor', rgb(col{ROI_exist(i)}),...
        'EdgeColor', [.2 .2 .2], 'LineWidth', .1);
end
set(gca, 'XTick', 1:length(uNames), 'XTickLabel', uNames,...
    'XLim', [0 length(uNames)+1], 'XTickLabelRotation', 90,...
    'FontSize', 5);
legend(Coord(1,:), 'Location', 'NorthEastOutside', 'Color', [.6 .6 .6]);
ylabel('relative frequency');
drawnow;

% plot reads
figure;
plotall(name(In), pos(In,:), background_image, scale);
drawnow;

% write output file
write_roi_countfile(fullfile(output_directory, 'DensityROICounts.csv'),...
    Coord, uNames, ROI_count, ROI_area)


% clustering (beta test version, only compatible with >=R2014b)
if length(ROI_exist) >= 3
    L = linkage(ROI_freq(idx_sort,:)','average');
    figure;
    sh1 = subplot('Position',[.05 .15 .05 .8]);
    dendrogram(L,'Orientation','left')
    set(sh1,'color','none',...
        'XLim',[-.001 inf],'YLim',[.5 length(ROI_exist)+.5]);
    idx_cluster = get(sh1,'YTickLabel');
    axis off
    idx_cluster = str2num(idx_cluster);
    idx_cluster = flipud(idx_cluster);
    
    sh2 = subplot('Position',[.1 .15 .8 .8]);
    bh = bar3(1:length(ROI_exist),ROI_freq(idx_sort,idx_cluster)',1);
    
    for i = 1:length(bh)
        Zdata = get(bh(i),'ZData');
        set(bh(i),'CData',Zdata);
    end
    YLab = Coord(1,:);
    YLab = YLab(idx_cluster);
    view(2)
    set(sh2,'color','none',...
        'XLim',[.5 length(idx_rest)+.5],'YLim',[.5 length(ROI_exist)+.5],...
        'plotboxaspectratiomode','auto',...
        'YAxisLocation','Right',...
        'YTick',1:length(ROI_exist),'YTickLabel',YLab,...
        'XTick',1:length(idx_rest),'XTickLabel',uNames(idx_sort),...
        'XTickLabelRotation',90, 'FontSize', 6);
end

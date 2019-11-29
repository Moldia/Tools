% draw polygons to create ROI 
% Polygons can be adjusted after they are created. When everything is done,
% double click to confirm the polygons
% Xiaoyan, 2018

clear;
close all;
drawnow;

%% modify here
decoded_file = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\3\QT_0.65_details_noNNNN.csv';
background_image = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\3\Sample 3_Base 3_resize20_c1.jpg';     % important for size
scale = .2;  % image scale

roi_number = 5;   % how many ROI regions you want to extract
col = {'white' 'yellow' 'blue' 'tomato' 'palevioletred' 'bisque' 'pink' 'mediumspringgreen'};

output_directory = 'E:\PROOOJECTS\test_dataset\ROI_draw_onImage';

%% do not modify

% load
[name, pos] = getinsitudata(decoded_file);
[name, pos] = removereads(name, 'NNNN', pos);
[uNames, ~, iName] = unique(name);
cName = hist(iName, 1:numel(uNames));

disp('loading image..')
image = imread(background_image);

posScaled = correctcoord(pos, scale);

% create ROIs and get polygon coordinates (in original scale)
[Coord, Coord_write, blank] = drawpolygon(image, roi_number, scale);

% remove empty polygons
if ~isempty(blank)
    Coord(:,blank) = [];
end

% write coordinates to file
mkdir(output_directory);
fid = fopen(fullfile(output_directory, 'ImageROICoordinates.csv'), 'w');
fprintf(fid,'Polygon id,x coordiates,y coordinates\n');
fprintf(fid,'%d,%d,%d\n', Coord_write');
fclose(fid);

% count transcripts within polygons
ROI_exist = 1:roi_number;
ROI_exist(blank) = [];
[ROI_count, ROI_freq, ROI_proportion, ROI_area, In] = ...
    count_in_polygon(length(ROI_exist), Coord,...
    pos, uNames, iName, cName);

% plot reads
plotall(name(In), pos(In,:), background_image, scale);
drawnow;

% plot polygons
if length(col)<roi_number
    col = repmat(col, 1, ceil(roi_number/length(col)));
end
figure;
imshow(image,[]);
hold on;
for i = 1:length(ROI_exist)
    plot(Coord{2,i}*scale, Coord{3,i}*scale,...
        'LineWidth', 2, 'Color', rgb(col{ROI_exist(i)}));
end
legend(Coord(1,:), 'Location', 'NorthEastOutside', 'Color', [.6 .6 .6]);
axis off;
drawnow;

% bar plot
figure;
bh = bar(ROI_freq);
set(gca, 'XTick', 1:numel(uNames), 'XTickLabel', uNames,...
    'XLim', [0 numel(uNames)+1], 'XTickLabelRotation', 90,...
    'FontSize', 6);
legend(Coord(1,:), 'Location', 'NorthEastOutside', 'color', [.6 .6 .6]);
ylabel('relative frequency');
box off;
for i = 1:length(bh)
    set(bh(i), 'FaceColor', rgb(col{ROI_exist(i)}), 'EdgeColor',[.2 .2 .2], 'LineWidth',0.1);
end
drawnow;

% write output file
write_roi_countfile(fullfile(output_directory, 'ImageROICounts.csv'),...
    Coord, uNames, ROI_count, ROI_area)

%% clustering (beta test version, only compatible with >=R2014b)
if length(ROI_exist) >= 3
    L = linkage(ROI_freq', 'average');
    figure;
    sh1 = subplot('Position', [.05 .15 .05 .8]);
    dendrogram(L, 'Orientation', 'left')
    set(sh1, 'color', 'none',...
        'XLim', [-.001 inf], 'YLim', [.5 length(ROI_exist)+.5]);
    idxcluster = get(sh1, 'YTickLabel');
    axis off
    idxcluster = str2num(idxcluster);
    idxcluster = flipud(idxcluster);
    
    L2 = linkage(ROI_freq, 'average');
    idxname = dendroperm(L2, numel(uNames));
    idxname = fliplr(idxname);
    
    sh2 = subplot('Position', [.1 .15 .8 .8]);
    bh = bar3(ROI_freq(idxname, idxcluster)', 1);
    for i = 1:length(bh)
        bh(i).CData = bh(i).ZData;
    end
    YLab = Coord(1,:);
    YLab = YLab(idxcluster);
    view(2)
    set(sh2, 'color', 'none',...
        'XLim', [.5 numel(uNames)+.5], 'YLim', [.5 numel(idxcluster)+.5],...
        'YAxisLocation', 'Right',...
        'YTick', 1:numel(ROI_exist), 'YTickLabel', YLab,...
        'XTick', 1:numel(uNames), 'XTickLabel', uNames(idxname),...
        'XTickLabelRotation', 90, 'FontSize', 6);
    axis normal
    grid on
end

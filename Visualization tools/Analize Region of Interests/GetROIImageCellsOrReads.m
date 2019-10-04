% create image subset containing ROI and extract reads or cells within polygon
% Xiaoyan, 2018

clear;
close all;

%% modify here
polygon = csvread('E:\PROOOJECTS\test_dataset\ROI_DrawDensity\DensityROICoordinates.csv', 1);
background_image = 'E:\PROOOJECTS\test_dataset\860502_1_align.png';
scale = .2;     % image scale
cells_or_reads_file = 'E:\PROOOJECTS\test_dataset\QT_0.35_0.004_details.csv';
output_directory = 'E:\PROOOJECTS\test_dataset\inROI';

%% do not modify

% load
data = importdata(cells_or_reads_file);
textcolumns = size(data.textdata,2) - size(data.data,2);

% find columns with positional data
poscolumns = cellfun(@(v) contains(lower(v), {'centroid', 'pos', 'location','global'}) & ~contains(lower(v), {'tilepos'}),...
    data.textdata(1,:));
if nnz(poscolumns)==2
    pos = data.data(:,find(poscolumns)-textcolumns);
else
    error('unexpected number of columns for position data');
end

% file names
[~, imname] = fileparts(background_image);
[~, fname] = fileparts(cells_or_reads_file);

roundedCoord = round(pos);
Iback = imread(background_image);

mkdir(output_directory);
for i = 1:max(polygon(:,1))
    coordinates = polygon(polygon(:,1)==i,2:3);
    
    if size(coordinates,1) > 3
        % a bit of buffer, can be excluded, used to create an image subset
        % containing only ROI+-5 px
        imin = max(floor(min(coordinates(:,2))-5), 1/scale);
        imax = min(ceil(max(coordinates(:,2))+5), floor(max(pos(:,2))));
        jmin = max(floor(min(coordinates(:,1))-5), 1/scale);
        jmax = min(ceil(max(coordinates(:,1))+5), floor(max(pos(:,1))));

        mask = poly2mask(coordinates(:,1)-jmin, coordinates(:,2)-imin, imax-imin+1, jmax-jmin+1);

        % square area
        insquare = readsinsqr(pos, [jmin imin jmax imax]);
        if ~nnz(insquare)
            break;
        end

        subpos = bsxfun(@minus, pos(insquare,:), [jmin, imin]);
        insquare = find(insquare);

        % polygon
        subpos(subpos<.5) = .5;
        inroi = logical(readsinroi(subpos, mask));
        inroi = insquare(inroi);

        clf; imshow(mask);
        hold on; plot(pos(inroi,1)-jmin, pos(inroi,2)-imin, '.');
        title(['ROI' num2str(i)]);
        drawnow;

        % save subset ROI images
        mask = imresize(mask, scale);
        Dapi = Iback(round(imin*scale):round(imax*scale),...
            round(jmin*scale):round(jmax*scale),:);
        
        mask = padimg(cast(mask, class(Iback)), size(Dapi,2)-size(mask,2), size(Dapi,1)-size(mask,1));
        % if color image, make mask 3D
        if length(size(Iback)) > 2
            mask = repmat(mask, 1, 1, 3);
        end
        
        Dapi = Dapi.*mask;
        clf; imshow(Dapi); title(['ROI' num2str(i)]); drawnow;
        imwrite(Dapi, fullfile(output_directory, [imname '_ROI' num2str(i) '.tif']));

        % write gene name and coordinate of reads in polygon, and original
        % index (coordinate system matching subset ROI image)
        % format kept same as input
        % one file for each ROI
        header = data.textdata(1,:);
        header = [header(~poscolumns), {'sub_x', 'sub_y', 'original_index'}];

        rownames = data.textdata(2:end,:);

        if ~isempty(inroi)
            if textcolumns
                towrite = [rownames(inroi,1:textcolumns), num2cell(data.data(inroi,~poscolumns(textcolumns+1:end))), num2cell([pos(inroi,1)-jmin, pos(inroi,2)-imin, inroi])]';
            else
                towrite = [num2cell(data.data(inroi,~poscolumns)), num2cell([pos(inroi,1)-jmin, pos(inroi,2)-imin, inroi])]';
            end
            
            fid = fopen(fullfile(output_directory, [fname '_ROI' num2str(i) '.csv']), 'w');
            fprintf(fid, lineformat('%s', numel(header)), header{:});
            fprintf(fid, [repmat('%s,', 1, textcolumns), lineformat('%d', size(data.data,2)+1)], towrite{:});
            fclose(fid);

        end
    end
end

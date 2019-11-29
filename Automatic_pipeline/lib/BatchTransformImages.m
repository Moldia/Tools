% batch transform images based on the given rotation and translation parameters
% Xiaoyan, 2017

%% input
transformation_file = 'registration.csv';
% columns: filename_prefix, 
reference_image_file = '';

%%
% reference image size (needed to make all image coordinates start from
% same origo)
szRef = [];
try 
    szRef = imfinfo(reference_image_file);
    szRef = [szRef.Height, szRef.Width];
end

% import registration parameters
reg_info = importdata(transformation_file);
reg_info = cell2table([reg_info.textdata(2:end,1:3), num2cell(reg_info.data)],...
    'VariableNames', strrep(reg_info.textdata(1,:), ' ', '_'));

for i = 1:size(reg_info, 1)     % image sets
    
    % find images to transform
    files = ls(reg_info.path{i});
    files = cellstr(files);
    files = files(cellfun(@(v)...
        strcmp(...
        v(1:min(length(v), length(reg_info.filename_prefix{i}))),...
        reg_info.filename_prefix{i}),...
        files));
    
    for j = 1:numel(files)      % images in a set
        I = imread(fullfile(reg_info.path{i}, files{j}));
        szFloat = size(I);
        
        name = strsplit(files{j}, filesep);
        name = strsplit(name{end}, '.tif');
        
        mkdir(fullfile(reg_info.path{i}, 'transformed'));

        name = fullfile(fullfile(reg_info.path{i}, 'transformed'),...
            [name{1}, '_transformed.tif']);
        
        try
            deltaX = reg_info.after_x(i) - reg_info.before_x(i) + ...
                (szRef(2) - szFloat(2))/2;  % anchor to upper left corner
            deltaY = reg_info.after_y(i) - reg_info.before_y(i) + ...
                (szRef(1) - szFloat(1))/2;  % anchor to upper left corner
            % image transform
            transformimage(I, -reg_info.rotation(i), ...
                -round(deltaX), -round(deltaY),...
                1, name, szRef);
        catch
            deltaX = reg_info.after_x(i) - reg_info.before_x(i);
            deltaY = reg_info.after_y(i) - reg_info.before_y(i);
            % image transform
            transformimage(I, -reg_info.rotation(i), ...
                -round(deltaX), -round(deltaY),...
                1, name);
        end
        
    end
end


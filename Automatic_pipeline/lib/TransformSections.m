% tranform reference image and spot coordinates based on given rotation and
% transformation parameters
% Xiaoyan, 2018

clear;
close all;

%% input
image_to_transform = 'Exp_2805_5_50_aligned.tif';
image_scale = .5;   % compared to coordinates
coordinates_file_to_transform = 'Exp_2805_5_QT_0.25_0.003_details.csv';
reference_image = 'Mut_2805_5_50_alignExp.tif';

% transformation parameters
scaling = 1;
rotation = 0;     % in degree, positive: counterclockwise
translation = [7860-7814, 5535-5431];     % shift in x and y, end - start
output_folder = 'Transformed';

%%
% organize data
I = imread(image_to_transform);
I = imresize(I, scaling);
imsize = size(I);

% transformation parameters
szRef = imfinfo(reference_image);
szRef = [szRef.Height, szRef.Width];
szFloat = imsize(1:2);
deltaX = translation(1) + (szRef(2) - szFloat(2))/2;
deltaY = translation(2) + (szRef(1) - szFloat(1))/2;

% transform image
mkdir(output_folder);
imfname = strsplit(image_to_transform, filesep);
imfname = strsplit(imfname{end}, '.tif');
imfname = fullfile(output_folder, [imfname{1}, '_Transformed.tif']);

if numel(imsize) > 2
    I2 = zeros(szRef(1), szRef(2), 3, class(I));
    for j = 1:3
        Itemp = transformimage(I(:,:,j),...
            rotation, -round(deltaY),  -round(deltaX), 1, imfname, szRef, 0);
        I2(:,:,j) = Itemp;
    end
    imwrite(I2, imfname);
else
    I2 = transformimage(I, rotation, -round(deltaY), -round(deltaX), 1, imfname, szRef);
end

% transform coordinates
[name, pos] = getinsitudata(coordinates_file_to_transform);
pos2 = transform_coordinates(pos, scaling, szFloat, rotation,...
     -translation(2), -translation(1), szRef, image_scale);

% save file
fname = strsplit(coordinates_file_to_transform, filesep);
fname = strsplit(fname{end}, '.csv');
fname = fullfile(output_folder, [fname{1}, '_Transformed.csv']);
replace_file_columns(coordinates_file_to_transform, [3,4], [{pos(:,1)},{pos(:,2)}], fname);

% visualization
figure; 
plotall(name, pos2, imfname, image_scale);




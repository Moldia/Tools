%% align a floating image on top of a reference image
%  need to give rotation angle and translation matrix
%  last update, 2015-2-18, Xiaoyan

clear;
close all;

%% input
ref_image = 'Stitched\base1_c1_stitched.tif'; % give sample snapshot image (blue DAPI)

input_image_prefix = 'Stitched\base3_c';
flo_image = [input_image_prefix '1_stitched.tif']; % give sample snapshot image (blue DAPI)

output_image_prefix = 'base3';

%% original
ref = imread(ref_image);
size_ref = size(ref);
flo = imread(flo_image);
ref_resized = imresize(ref, .1);
clear ref;
flo_resized = imresize(flo, .1);
Ifuse = imfuse(flo_resized, ref_resized);
imshow(Ifuse);
% green: floating
% purple: reference

%% rotation
angle = 0; % positive: counter clockwise
[flo_rotate, Ifuse_rotate] = rotateimage(flo_resized, angle, ref_resized);
imshow(Ifuse_rotate);

%% translation
yup = -2;   % positive: move the floating image up, negative: down
xleft = 8;   % positive: move the floating image left, negative: right
Ifuse_translate = translateimage(yup, xleft, flo_rotate, ref_resized);
imshow(Ifuse_translate);

%% transform images
mkdir('AlignedImages');
transformimage(flo,angle,yup,xleft,10,...
    ['AlignedImages\' output_image_prefix '1_ORG.tif'],size_ref);
tic
for c = 2:6
    flo = imread([input_image_prefix num2str(c) '_ORG.tif']);
    transformimage(flo,angle,yup,xleft,10,...
        ['AlignedImages\' output_image_prefix num2str(c) '_ORG.tif'],size_ref);
end
toc
%% align a floating image on top of a reference image
%  need to give rotation angle and translation matrix
%  last update, 2015-2-18, Xiaoyan

clear;
close all;

%% input
ref_image = 'F:\191001_OL\191004_mBrain_combined\1\Aligned_Images_Rigid\Test\1\Aligned_Images_Rigid\Base 1_c1_ORG.tif'; % give sample snapshot image (blue DAPI)

input_image_prefix = 'F:\191001_OL\191004_mBrain_combined\1\Aligned_Images_Rigid\Test\1\Aligned_Images_Rigid\Base 2_c';
flo_image = [input_image_prefix '1_ORG.tif']; % give sample snapshot image (blue DAPI)

output_image_prefix = 'F:\191001_OL\191004_mBrain_combined\1\Aligned_Images_Rigid\Test\1\Aligned_Images_Rigid\Base 3_c';

%% original
ref = imread(ref_image);
size_ref = size(ref);
flo = imread(flo_image);
ref_resized = imresize(ref, .5);
clear ref;
flo_resized = imresize(flo, .5);
Ifuse = imfuse(flo_resized, ref_resized);
imshow(Ifuse);
% green: floating
% purple: reference

%% rotation
angle = 0; % positive: counter clockwise
[flo_rotate, Ifuse_rotate] = rotateimage(flo_resized, angle, ref_resized);
imshow(Ifuse_rotate);

%% translation
yup = 0;   % positive: move the floating image up, negative: down
xleft = 0;   % positive: move the floating image left, negative: right
Ifuse_translate = translateimage(yup, xleft, flo_rotate, ref_resized);
imshow(Ifuse_translate);

%% transform images
%mkdir('D:\170717_schizo_CA1Probes\exported\1441-hippo\aligned\170717_mBrain_schizo_1441-hippo_b1');
transformimage(flo,angle,yup,xleft,2,...
    [output_image_prefix '1_ORG.tif'],size_ref);
tic
for c = 2:6
    flo = imread([input_image_prefix num2str(c) '_ORG.tif']);
    %mkdir([output_image_prefix num2str(c)]);
    transformimage(flo,angle,yup,xleft,2,...
        [output_image_prefix num2str(c) '_ORG.tif'],size_ref);
end
toc
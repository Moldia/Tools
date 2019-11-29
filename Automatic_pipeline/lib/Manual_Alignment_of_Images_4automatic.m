%% align a floating image on top of a reference image
%  need to give rotation angle and translation matrix
%  last update, 2015-2-18, Xiaoyan


%% input
ref_image = 'E:\Melanoma\Align\Base_1_aligned-1.tif'; % give sample snapshot image (blue DAPI)

input_image_prefix = 'E:\Melanoma\Align\Base_4_aligned-';
flo_image = [input_image_prefix '1.tif']; % give sample snapshot image (blue DAPI)

output_image_prefix = flo_image;

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
yup = -2;   % positive: move the floating image up, negative: down
xleft = 0;   % positive: move the floating image left, negative: right
Ifuse_translate = translateimage(yup, xleft, flo_rotate, ref_resized);
imshow(Ifuse_translate);

%% transform images
%mkdir('D:\170717_schizo_CA1Probes\exported\1441-hippo\aligned\170717_mBrain_schizo_1441-hippo_b1');
transformimage(flo,angle,yup,xleft,2,...
    [output_image_prefix '1.tif'],size_ref);
tic
for c = 2:6
    flo = imread([input_image_prefix num2str(c) '.tif']);
    %mkdir([output_image_prefix num2str(c)]);
    transformimage(flo,angle,yup,xleft,2,...
        [output_image_prefix num2str(c) '.tif'],size_ref);
end
toc
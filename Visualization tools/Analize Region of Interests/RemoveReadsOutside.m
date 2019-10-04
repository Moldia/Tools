% remove reads outside tissue
% Xiaoyan, 2017


I = imread('slideA\DO_c1.tif');
scale = 1;

I = imfilter(I, fspecial('average', 10));
Ibw = im2bw(I, graythresh(I));
Ibw = imdilate(Ibw, strel('disk',3));
Ibw = imfill(Ibw, 'holes');

% visualize segmented mask and wait until the figure window is closed
Ibound = imresize(Ibw,.2) ~= imerode(imresize(Ibw,.2), strel('disk',1));
Icolor = repmat(imresize(uint8(double(I)/double(max(I(:)))*225), .2), 1, 1, 3);
Icolor(:,:,1) = Icolor(:,:,1) + uint8(Ibound)*255;
h = imshow(Icolor);
% h = imshow(Ibw, []);
waitfor(h);

% save mask
imwrite(Ibw, ['TissueMask_' num2str(scale*100) '%.tif']);

% save coordinates
fname = savereadsinbw('Decoding\QT_0.45_0_details.csv', Ibw, scale);

% plot
[name, pos] = getinsitudata(fname);
[name, pos] = removereads(name, 'NNNN', pos);
plotall(name, pos, ['TissueMask_' num2str(scale*100) '%.tif'], scale);

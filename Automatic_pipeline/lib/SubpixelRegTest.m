% comparison of two image registration methods
% imregcorr (MATLAB image processing toolbox) and DFTregistration
% Xiaoyan, 2017

%% 
I1 = imread('b1_c2.tif');
I2 = imread('b3_c2.tif');

Ax = [];

ax = subplot(221); 
Ax = [Ax, ax];
imshowpair(I1, I2);
axis on
title('original');


%% imregcorr
tic
tform = imregcorr(I1, I2);
tform.T

% pixel-resolution translation only
I2reg = padimg(I2, -round(tform.T(3,1)), -round(tform.T(3,2)), 'NW');
toc

ax = subplot(222); 
Ax = [Ax, ax];
imshowpair(I1, I2reg);
axis on
title('imregcorr');

%% DFT with pixel resolution
tic 
[tform, Greg] = dftregistration(fft2(I1), fft2(I2), 1);
tform

I2reg = padimg(I2, round(tform(4)), round(tform(3)), 'NW');
toc

ax = subplot(223); 
Ax = [Ax, ax];
imshowpair(I1, I2reg);
axis on
title('dftregistration, no upscaling');

% inverse FFT
tic
I2regIFFT = ifft2(Greg);
I2regIFFT = uint16(real(I2regIFFT));
toc

%% DFT registration with 0.05 subpixel resolution
tic 
[tform, Greg] = dftregistration(fft2(I1), fft2(I2), 20);
tform

I2reg = padimg(I2, round(tform(4)), round(tform(3)), 'NW');
toc

% inverse FFT
tic
I2regIFFT = ifft2(Greg);
I2regIFFT = uint16(real(I2regIFFT));
toc

ax = subplot(224); 
Ax = [Ax, ax];
imshowpair(I1, I2regIFFT);
axis on
title('dftregistration, 20x upscaling');
 
%%
linkaxes(Ax, 'xy');




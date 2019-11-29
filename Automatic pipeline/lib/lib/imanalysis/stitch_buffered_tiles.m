function I = stitch_buffered_tiles...
    (refimfile, tileprefix, centersize, buffersize, tilesuffix, ndigits, imgnum)
% I = stitch_buffered_tiles(refimfile, tileprefix, centersize, buffersize, tilesuffix, ndigits, imgnum)
% take only center part of aligned tiles and stitch
% Xiaoyan, 2017


imsize = imfinfo(refimfile);
imsize = [imsize.Height, imsize.Width];

ntileX = ceil(imsize(2)/centersize);
ntileY = ceil(imsize(1)/centersize);

if nargin <= 4
    tilesuffix = '.tif';
end

isIcreated = 0;

if nargin < 5
    ndigits = ceil(log10(ntileX*ntileY));
end

for i = 1:ntileY
    for j = 1:ntileX
        try
            tnum = paddigits(imgnum(ntileX*(i-1)+j), ndigits);
        catch
            tnum = paddigits((ntileX*(i-1)+j), ndigits);
        end
        im = imread([tileprefix, tnum, tilesuffix]);
        if ~isIcreated
            try
                I = zeros(ntileY*centersize, ntileX*centersize, size(im,3), class(im));
            catch
                I = zeros(ntileY*centersize, ntileX*centersize, 1, class(im));
            end
            isIcreated = 1;
        end
        imcrop = im(buffersize+1:buffersize+centersize,...
            buffersize+1:buffersize+centersize, :);
        I(centersize*(i-1)+1:centersize*i,...
            centersize*(j-1)+1:centersize*j, :) = imcrop;
    end
end

I = padimg(I, imsize(2)-size(I,2), imsize(1)-size(I,1));

end
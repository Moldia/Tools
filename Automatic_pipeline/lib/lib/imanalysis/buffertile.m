function buffertile(im, centersize, buffersize, outdir)
% buffertile(im, centersize, buffersize, outdir)
% tiles with buffer
% Xiaoyan, 2017


% pad image
imsize = size(im);
ntileX = ceil(imsize(2)/centersize);
ntileY = ceil(imsize(1)/centersize);
im = padimg(im,...
    centersize*ntileX-imsize(2)+buffersize,...
    centersize*ntileY-imsize(1)+buffersize, 'SE');
im = padimg(im, buffersize, buffersize, 'WN');

% tile and save
mkdir(outdir);
for i = 1:ntileY
    for j = 1:ntileX
        xstart = centersize*(j-1);
        ystart = centersize*(i-1);
        tile = im(ystart+1:ystart+centersize+2*buffersize,...
            xstart+1:xstart+centersize+2*buffersize,:);
        imwrite(tile,...
            [outdir, '\tile' num2str(ntileX*(i-1)+j), '.tif']);
    end
end

end

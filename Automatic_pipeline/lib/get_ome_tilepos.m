function [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] =...
    get_ome_tilepos(reader)
% [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] =...
%     get_ome_tilepos(reader)
% input: OMERO reader or microscopy image file that contains OMERO metadata
% Xiaoyan, 2017
% get reader if input is file name (character)
if ischar(reader)
    reader = bfGetReader(reader);
end

nSeries = reader.getSeriesCount();
nImage = reader.getImageCount();
nChannels = reader.getSizeC;
omeMeta = reader.getMetadataStore();

%disp(reader.getSeriesMetadataValue());

nZstacks = omeMeta.getPixelsSizeZ(0).getValue();
pixelsize = double(omeMeta.getPixelsPhysicalSizeX(0).value());

xypos  = [];
for i = 1:nSeries
    if isempty(omeMeta.getPlanePositionX(i-1, 0)) || isempty(omeMeta.getPlanePositionY(i-1, 0))
        %fprintf('Could not load tile positions for series %d\n', i);
    else
        posX = double(omeMeta.getPlanePositionX(i-1, 0).value());
%         disp(posX);
%         pause();
        posY = double(omeMeta.getPlanePositionY(i-1, 0).value());
        xypos = [xypos; posX, posY];
        fprintf('Loaded tile positions for series %d\n', i);
    end
end

if ~isempty(xypos)
    xypos = bsxfun(@minus, xypos, min(xypos,[],1))/pixelsize;
end

nSerieswPos = size(xypos,1);
%nSerieswPos=reader.getImageCount()
% close
if ischar(reader)
    reader.close();
end

end

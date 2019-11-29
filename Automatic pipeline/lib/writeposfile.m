function writeposfile(refimf, tilesize, hyb, channels, refchannels,...
    hybfolderprefix, channelfolderprefix, channelfoldersuffix, outputfolder)
% writeposfile(refimf, tilesize, hyb, channels, refchannels,...
%     hybfolderprefix, channelfolderprefix, channelfoldersuffix, outputfolder)
% simply write position file for CP input without processing images
% Xiaoyan, 2017


imsize = imfinfo(refimf);
imsize = [imsize.Width, imsize.Height];
ntiles = ceil(imsize/tilesize);

% get image directory names
imdirs = cell(ntiles(1)*ntiles(2), length(hyb), length(channels));
for b = 1:length(hyb)
    hybnum = hyb(b);
    for c = 1:length(channels)
        cnum = channels(c);
        [imdir, imfiles] = imfiles_in_grid(...
            [hybfolderprefix, num2str(hybnum),...
            channelfolderprefix, num2str(cnum), channelfoldersuffix],...
            'tile', '.tif', 0, ntiles);
        imdirs(:,b,c) = reshape(imdir', [], 1);
    end
end
imfiles = reshape(imfiles', [], 1);

% tile starting position
tileposX = repmat(0:tilesize:tilesize*(ntiles(1)-1), ntiles(2), 1);
tileposY = repmat((0:tilesize:tilesize*(ntiles(2)-1))', 1, ntiles(1));

% per hyb data
towrite = [];
header = catstrnum('channel', channels);
metadata = num2cell([(1:ntiles(1)*ntiles(2))', reshape(tileposX',[],1), reshape(tileposY',[],1)]);
for b = 1:length(hyb)
    perbase = [];
    for c = 1:length(channels)
        perbase = [perbase, imdirs(:,b,c), imfiles(:)];
    end
    towrite = [towrite; [metadata, num2cell(repmat(b, ntiles(1)*ntiles(2), 1)), perbase]];
end

% add reference channels (default 1st in hyb)
for i = 1:length(refchannels)
    towrite = [towrite,...
        repmat([imdirs(:,1,channels == refchannels(i)), imfiles], length(hyb), 1)];
    header = [header, ['General_channel' num2str(refchannels(i))]];
end
[~, idx] = sort(cell2mat(towrite(:,1)));
towrite = towrite(idx,:)';

% file header
header = [strcat('Image_PathName_', header); strcat('Image_FileName_', header)];
header = [{'Metadata_position','Tile_xPos','Tile_yPos','Hyb_step'}, reshape(header, 1, [])];

% write file
mkdir(outputfolder)
fid = fopen(fullfile(outputfolder, 'Tiled.csv'), 'w');
fprintf(fid, lineformat('%s', length(header)), header{:});
fprintf(fid, ['%d,%d,%d,%d,' lineformat('%s', length(header)-4)], towrite{:});
fclose(fid);

end


    
% extended depth of focus, then tophat
% Xiaoyan, 2017

mkdir('Preprocessing');
nChannels = 6;
nSeries = 162;
nZstacks = 7;

for c = 1:nChannels
    if c == 1
        SE = strel('disk', 25);     % DAPI
    else
        SE = strel('disk', 3);
    end
    
    for s = 1:nSeries
        if mod(s, 10) == 0
            fprintf('Channel %d. Tiles finished: %d out of %d\n', c, s, nSeries);
        end
        tnum = paddigits(s, 3);
        
        I = cell(nZstacks,1);
        parfor z = 1:nZstacks
            im = imread(['ZENout\base3\base3_z', num2str(z),...
                'c', num2str(c), 'm', tnum, '_ORG.tif']);
            I{z} = im;
        end
        
        % focus stacking and tophat
        IFS = fstack_modified(I);
        IFS = imtophat(IFS, SE);
        
        % write stack image
        imwrite(IFS, ['Preprocessing\base3', ...
            '_t', num2str(s), '_FS_tophat_stack.tif'], 'tiff',...
            'writemode', 'append');
    end
end





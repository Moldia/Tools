% extended depth of focus, then tophat
% Xiaoyan, 2017


for s = 1:1     % series
    output_dirname = 'FS';
    output_prefix = 'base1';
    
    imfile = ['D:\200_12_2_base' num2st(s) '.czi'];
    outfolder = ['D:\200_12_2_base1\Preprocessing\', output_dirname];
    mkdir(outfolder);
    
    % Construct a Bio-Formats reader decorated with the Memoizer wrapper
    bfreader = loci.formats.Memoizer(bfGetReader(), 0);
    % Initialize the reader with an input file
    bfreader.setId(imfile);
    
    [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
        get_ome_tilepos(bfreader);
    scene = nSeries/nSerieswPos;
    
    bfreader.close();
    
    for c = 1:nChannels
        if c == 1
            SE = strel('disk', round(8/pixelsize));     % DAPI
        else
            SE = strel('disk', round(1/pixelsize));
        end
        
        parfor t = 1:nSerieswPos
            fprintf('Channel %d. Tile %d / %d\n',...
                c, t, nSerieswPos);
            
            % Initialize a new reader per worker as Bio-Formats is not thread safe
            bfreader2 = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
            % Initialization should use the memo file cached before
            bfreader2.setId(imfile);
            
            bfreader2.setSeries(scene*t-1);
            
            I = cell(nZstacks,1);
            for z = 1:nZstacks
                iPlane = bfreader2.getIndex(z-1, c-1, 0)+1;
                I{z} = bfGetPlane(bfreader2, iPlane);
            end
            bfreader2.close();

            % focus stacking
            IFS = fstack_modified(I);
            
%             % or MIP
%             IFS = max(cat(3, IFS{:}), [], 3);
            
            % tophat
            IFS = imtophat(IFS, SE);
            
            % write stack image
%             tnum = paddigits(t, 3);
            imwrite(IFS, fullfile(outfolder,...
                [output_prefix, '_t', num2str(t), '_', output_dirname, '.tif']),...
                'tiff', 'writemode', 'append');
        end
        datetime
    end
    
    % tile position configuration file (for stitching)
    fid = fopen(fullfile(outfolder, ['TileConfiguration_', output_prefix, '.txt']), 'w');
    fprintf(fid,'# Define the number of dimensions we are working on\n');
    fprintf(fid,'dim = 2\n\n# Define the image coordinates\n');
    for t = 1:nSerieswPos
        fprintf(fid,...
            '%s_t%d_%s.tif; ; (%.1f, %.1f)\n',...
            output_prefix, t, output_dirname, xypos(t,1), xypos(t,2));
    end
    fclose(fid);
    
    % or run MIST (separate script)
    
end


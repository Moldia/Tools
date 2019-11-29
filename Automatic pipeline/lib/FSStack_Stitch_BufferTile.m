% extended depth of focus, tophat and stitch, and tile
% Xiaoyan, 2017

foldernames = {'F:\161220KI', 'F:\161220KI', 'F:\160519KI'};
sectionnames = {'170315_161220KI_4-1', '170315_161220KI_4-3', '170315_160519KI_B_1'};
input_prefix = '_b1';

for s = 1:length(sectionnames)
    output_dirname = 'FS';
    output_prefix = 'base1';
    
    imfile = fullfile(foldernames{s}, [sectionnames{s}, input_prefix '.czi']);
    outfolder = fullfile(sectionnames{s}, 'Preprocessing', output_dirname);
    mkdir(outfolder);

    % Construct a Bio-Formats reader decorated with the Memoizer wrapper
    reader = loci.formats.Memoizer(bfGetReader(), 0);
    % Initialize the reader with an input file
    reader.setId(imfile);

    [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
        get_ome_tilepos(reader);
    scene = nSeries/nSerieswPos;
    
    reader.close();

    for c = 1:nChannels
%         if c == 1
%             SE = strel('disk', round(8/pixelsize));     % DAPI
%         else
%             SE = strel('disk', round(1/pixelsize));
%         end
        
        parfor t = 1:nSerieswPos
            fprintf('Channel %d. Tile %d / %d\n',...
                c, t, nSerieswPos);
            
            % Initialize a new reader per worker as Bio-Formats is not thread safe
            reader2 = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
            % Initialization should use the memo file cached before
            reader2.setId(imfile);
            
            reader2.setSeries(scene*t-1);
            
            I = cell(nZstacks,1);
            for z = 1:nZstacks
                iPlane = reader2.getIndex(z-1, c-1, 0)+1;
                I{z} = bfGetPlane(reader2, iPlane);
            end
            reader2.close();

            % focus stacking
            IFS = fstack_modified(I);
            
%             % tophat
%             IFS = imtophat(IFS, SE);

            % write stack image
%             tnum = paddigits(t, 3);
            imwrite(IFS,...
                fullfile(outfolder,...
                [output_prefix, ...
                '_t', num2str(t), '_', output_dirname, '.tif']),...]\
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
    
    % MIST stitching
    imwrite(zeros(2048, 'uint16'), fullfile(outfolder, 'empty.tif'));

    [~, img_name_grid] = tilepos2grid(xypos, [output_prefix '_t'], ['_' output_dirname '.tif'], 3);
    input_directory = outfolder;
    output_directory = fullfile(outfolder, 'Stitched');
    save_stitched_image = 1;
    assemble_from_metadata = 1;

    mkdir(output_directory);

    for c = 1:6
        if c == 1
            assemble_from_metadata = 0;
        else
            assemble_from_metadata = 1;
        end
        stitch_time_slice(input_directory, img_name_grid,...
            output_directory, output_prefix,...
            1, NaN, NaN,...     % these will probably not change
            'Max', 0,...        % blending method
            save_stitched_image, assemble_from_metadata,...
            fullfile(output_directory, ['log_', output_prefix '.txt']),...
            10, 10, c);
        movefile(fullfile(output_directory, [output_prefix, 'stitched-1.tif']),...
            fullfile(output_directory, [output_prefix, '_c' num2str(c) '_stitched.tif']));
    end
    
    % tiling with some buffer
    for c = 1:6
        I = imread(fullfile(output_directory, [output_prefix, '_c' num2str(c) '_stitched.tif']));
        buffertile(I,...
            2000, 100,...
            fullfile(output_directory, ['BufferTiled_', output_prefix, '_c' num2str(c)]));
    end
end


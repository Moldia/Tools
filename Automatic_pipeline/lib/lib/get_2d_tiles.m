%% get MIPs or make EDFs from z stacks
function [outdir, nFiles] = get_2d_tiles(indir, files, outdir)
disp('Starting get_2d_tiles')
outdir = fullfile(outdir, '2DTiles');
mkdir(outdir);

nFiles = zeros(numel(files),1);
disp(nFiles);
pause();
for s = 1:numel(files)
    s=1
    imfile = files{s};
    [~, output_prefix] = fileparts(imfile);
    
    % Construct a Bio-Formats reader decorated with the Memoizer wrapper
    bfreader = loci.formats.Memoizer(bfGetReader(), 0);
    % Initialize the reader with an input file
    bfreader.setId(fullfile(indir, imfile));
    
    [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
        get_ome_tilepos(bfreader);
    scene = nSeries/nSerieswPos;
    
    % check first available tile
    bfreader.setSeries(scene-1);
    tilesize = [bfreader.getSizeY, bfreader.getSizeX];
    bfreader.close();
    
    if any(tilesize ~= 2048)
        % tile size in image file does not match known FOV size
        warning('%s Tile size incorrect: X %d Y %d. Skipping..\n', imfile,  tilesize(2), tilesize(1));
    else
        nFiles(s) = nSerieswPos;
        parfor t = 1:nSerieswPos
            
            % skip if already exists
            if ~exist(fullfile(outdir, [output_prefix, '_t', num2str(t) '.tif']), 'file')
                
                % Initialize a new reader per worker as Bio-Formats is not thread safe
                bfreader2 = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
                % Initialization should use the memo file cached before
                bfreader2.setId(fullfile(indir, imfile));
                
                bfreader2.setSeries(scene*t-1);
                
                
                fprintf('%s Tile %d / %d\n', imfile, t, nSerieswPos);
                
                for c = 1:nChannels
                    if nZstacks == 1
                        % only one z plane in the file
                        iPlane = bfreader2.getIndex(0, c-1, 0)+1;
                        I = bfGetPlane(bfreader2, iPlane);
                        imwrite(I, fullfile(outdir,...
                            [output_prefix, '_t', num2str(t) '.tif']),...
                            'tiff', 'writemode', 'append');
                    else
                        I = cell(nZstacks,1);
                        for z = 1:nZstacks
                            iPlane = bfreader2.getIndex(z-1, c-1, 0)+1;
                            I{z} = bfGetPlane(bfreader2, iPlane);
                        end
                        % focus stacking
                        IFS = fstack_modified(I);
                        imwrite(IFS, fullfile(outdir,...
                            [output_prefix, '_t', num2str(t) '.tif']),...
                            'tiff', 'writemode', 'append');
                    end
                end
                bfreader2.close();
            end
            
        end
        datetime
    end
    
    % tile position configuration file
    csvwrite(fullfile(outdir, ['tile_coordinates_' output_prefix '.csv']), xypos);
    
    fid = fopen(fullfile(outdir, ['TileConfiguration_', output_prefix, '.txt']), 'w');
    fprintf(fid,'# Define the number of dimensions we are working on\n');
    fprintf(fid,'dim = 2\n\n# Define the image coordinates\n');
    for t = 1:nSerieswPos
        fprintf(fid,...
            '%s_t%d.tif; ; (%.1f, %.1f)\n',...
            output_prefix, t, xypos(t,1), xypos(t,2));
    end
    fclose(fid);
end
end
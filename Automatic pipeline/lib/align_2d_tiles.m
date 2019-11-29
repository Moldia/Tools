%% align all tiles to the reference round



function align_2d_tiles(indir, files, outdir, nChannels, nTiles, ref_round, ref_channel)
outdir = fullfile(outdir, 'Aligned2DTiles');
mkdir(outdir);

TForm = cell(nTiles, 1);
parfor t = 1:nTiles
    disp(strcat('Aligning tile',num2str(t)))
    Tform = zeros(numel(files),2);
    for s = [ref_round, setdiff(1:numel(files), ref_round)]
        imfile = files{s};
        [~, output_prefix] = fileparts(imfile);
        
        if s == ref_round
            % skip if already exists
            if ~exist(fullfile(outdir, [output_prefix, '_t', num2str(t) '.tif']), 'file')
                for c = 1:nChannels
                    I = imread(fullfile(indir, [output_prefix '_t' num2str(t) '.tif']), c);
                    imwrite(I, fullfile(outdir, [output_prefix '_t' num2str(t) '.tif']), 'WriteMode', 'append');
                    if c == ref_channel
                        refimg = I;
                        Tform(s,:) = [0 0];
                    end
                end
            else
                refimg = imread(fullfile(indir, [output_prefix '_t' num2str(t) '.tif']), ref_channel);
                Tform(s,:) = [0 0];
            end
			
        else
            for c = [ref_channel, setdiff(1:nChannels, ref_channel)]
                I = imread(fullfile(indir, [output_prefix '_t' num2str(t) '.tif']), c);
                if c == ref_channel
                    floatimg = I;
                    
                    % 0.5px-resolution registration
                    [tform, Greg] = dftregistration(fft2(refimg), fft2(floatimg), 2);
                    Tform(s,:) = [tform(4), tform(3)];
                end
                
                % but pixel-resolution transformation (no rotation)
                I = padimg(I, round(tform(4)), round(tform(3)), 'NW');
                I = padimg(I, 2048-size(I,2), 2048-size(I,1));
                imwrite(I, fullfile(outdir, [output_prefix '_t' num2str(t) '.tif']), 'WriteMode', 'append');
            end
        end
    end
    TForm{t} = Tform;
end
save(fullfile(outdir, 'transformation.mat'), 'TForm');
end

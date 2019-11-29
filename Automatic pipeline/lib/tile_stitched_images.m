%% re-slice images
function nTiles = tile_stitched_images(indir, files, outdir, tile_size)

outdir = fullfile(outdir, 'ReslicedTiles');
mkdir(outdir)

Stitched = ls(indir);
Stitched = cellstr(Stitched(3:end,:));

stitched_files = cell(0, numel(files));
for s = 1:numel(files)
    [~, filename] = fileparts(files{s});
    stitched_files{s} = Stitched(contains(Stitched, filename) & contains(Stitched, '_stitched'));
end
stitched_files = cat(2, stitched_files{:});
try
    imsize = cellfun(@(v) imfinfo(v), fullfile(indir, stitched_files(1,:)), 'UniformOutput', 0);
    imsize = cellfun(@(v) [v.Height, v.Width], imsize, 'UniformOutput', 0);
    imsize = cat(1, imsize{:});
catch ME
    if contains(ME.identifier, 'whatFormat')
        imsize = zeros(numel(files), 2);
    end
end

for s = 1:numel(files)
    for c = 1:size(stitched_files,1)
        try
            I = imread(fullfile(indir, stitched_files{c,s}));
        catch ME
            if contains(ME.identifier, 'fileFormat')
                load(fullfile(indir, stitched_files{c,s}));
            end
        end
                
        if ~isequal(ceil(min(imsize)/tile_size), ceil(imsize(s,:)/tile_size))
            I = padimg(I, min(imsize(:,2))-size(I,2), min(imsize(:,1))-size(I,1));
        end
        [~, output_prefix] = fileparts(stitched_files{c,s})
        output_prefix = strrep(output_prefix, '.tif', '');
        OUTNAME=['Base ' num2str(s) '_c' num2str(c) '_ORG']
        tileimage(I, tile_size, fullfile(outdir,OUTNAME));
    end
end

nX = ceil(size(I,2)/tile_size);
nY = ceil(size(I,1)/tile_size);
tileposX = repmat(0:tile_size:tile_size*(nX-1), nY, 1)';
tileposY = repmat((0:tile_size:tile_size*(nY-1))', 1, nX)';
nTiles = nX*nY;

csvwrite(fullfile(outdir, 'tilepos.csv'), [tileposX(:), tileposY(:)]);

end


image_directory = 'H:\Sample 1\';
image_files = {
    '190627_org1_cyc1(1).czi',
    '190629_org1_cyc2(1).czi',
    '190701_org1_cyc3(1).czi',
    '190702_org1_cyc4(1).czi',
    '190715_org1_cyc5(1).czi',
    
    };
output_directory = 'PreprocessingERIKOrg_v2';
nChannels = 6;
reference_round = 1;
tile_align_channel = 1; % for registration between rounds
tile_stitch_channel = 1; %Cy3   % for registration between tiles

% extract tile images from Zeiss files
[OUT, nFiles] = get_2d_tiles(image_directory, image_files, output_directory);

% align between rounds (assuming position offset is relatively small)
align_2d_tiles(fullfile(output_directory, '2DTiles'), image_files, output_directory, nChannels, max(nFiles), reference_round, tile_align_channel);

% stitch into big images either individually for every round or compute
% from a fixed round (need tile transformation info from align_2d_tiles)
stitch_2d_mist(fullfile(output_directory, '2DTiles'), image_files, output_directory, nChannels, tile_stitch_channel, reference_round);

% reslice into non-overlapping small tiles
%%nTiles = tile_stitched_images(fullfile(output_directory, 'Stitched2DTiles_MIST_Ref1'), image_files, output_directory, 2000);

% pick random 10 tiles to visually inspect (will pause until the figure
% window is closed)
%%inspect_images(fullfile(output_directory, 'ReslicedTiles'), image_files, nTiles, 2:5)

% preliminary base calling
% DO-1:4 corresponding to channel 4,5,2,3
%%Results = test_analysis(fullfile(output_directory, 'ReslicedTiles'), image_files, output_directory, nTiles, [4 5 2 3], 'taglist_mouse120.csv');



%% get MIPs or make EDFs from z stacks
function [outdir, nFiles] = get_2d_tiles(indir, files, outdir)
outdir = fullfile(outdir, '2DTiles');
mkdir(outdir);

nFiles = zeros(numel(files),1);

for s = 1:numel(files)
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




%% align all tiles to the reference round
function align_2d_tiles(indir, files, outdir, nChannels, nTiles, ref_round, ref_channel)
outdir = fullfile(outdir, 'Aligned2DTiles');
mkdir(outdir);

TForm = cell(nTiles, 1);
parfor t = 1:nTiles
    t
    
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



%% iss_suite stitching, modified from @iss/register
% results less optimal, not yet tested extensively
function stitch_2d_iss(indir, files, outdir, nChannels, ref_channel, corr_threshold, ref_file)

if nargin < 6
    corr_threshold = [.3 .6];
end

outdir_top = outdir;
if nargin < 7 || ref_file==0
    ref_file = 0;
    outdir = fullfile(outdir, 'Stitched2DTiles_ISS_Individual');
else
    outdir = fullfile(outdir, ['Stitched2DTiles_ISS_Ref' num2str(ref_file)]);
end

mkdir(outdir);

s_order = 1:numel(files);
if ref_file
    s_order = [ref_file, setdiff(s_order, ref_file)];
end

for s = s_order
    
    imfile = files{s};
    [~, output_prefix] = fileparts(imfile);
    
    xypos = csvread(fullfile(outdir_top, '2DTiles', ['tile_coordinates_' output_prefix '.csv']));
    
    
    [gridorder, tilefiles] = tilepos2grid(xypos,...
        [output_prefix '_t'], '.tif', 0);
    
    [nY, nX] = size(tilefiles);
    nTiles = nY*nX;
    NonemptyTiles = find(~strcmp(tilefiles, 'empty.tif'))';
    EmptyTiles = strcmp(tilefiles, 'empty.tif');
    
    if ~ref_file || s == ref_file
        
        try
            TileOrigin = csvread(fullfile(outdir, ['tile_origin_' output_prefix, '.csv']));
        catch
            
            RefImages = zeros(2048, 2048, nY, nX, 'uint16');
            
            for t = NonemptyTiles(:)'
                [y,x] = ind2sub([nY nX], t);
                if mod(t,10)==0; fprintf('Loading tile %d reference image\n', t); end
                Im = imread(fullfile(indir, tilefiles{y,x}), ref_channel);
                
                RefImages(:,:,t) = imfilter(Im, fspecial('disk', 1));
            end
            
            VerticalPairs = zeros(0,2);
            HorizontalPairs = zeros(0,2);
            vShifts = zeros(0,2);
            hShifts = zeros(0,2);
            ccv = zeros(0,1);
            cch = zeros(0,1);
            
            for t = NonemptyTiles
                [y,x] = ind2sub([nY nX], t);
                
                % align to south neighbor
                if y<nY && ~EmptyTiles(t+1)
                    [shift, cc, RefFFTStore] = ImRegFft2(RefImages(:,:,t), RefImages(:,:,t+1), corr_threshold, 200^2);
                    if all(isfinite(shift))
                        VerticalPairs = [VerticalPairs; t, t+1];
                        vShifts = [vShifts; shift];
                        ccv = [ccv; cc];
                    end
                    fprintf('%d, %d, down: shift %d %d, cc %f\n', y, x, shift, cc);
                end
                
                % align to east neighbor
                if x<nX && ~EmptyTiles(t+nY)
                    if isvarname(RefFFTStore)
                        [shift, cc] = ImRegFft2(RefImages(:,:,t), RefImages(:,:,t+nY), corr_threshold, 200^2);
                    else
                        [shift, cc] = ImRegFft2(RefImages(:,:,t), RefImages(:,:,t+nY), corr_threshold, 200^2);
                    end
                    if all(isfinite(shift))
                        HorizontalPairs = [HorizontalPairs; t, t+nY];
                        hShifts = [hShifts; shift];
                        cch = [cch; cc];
                    end
                    fprintf('%d, %d, right: shift %d %d, cc %f\n', y, x, shift, cc);
                    
                end
            end
            
            
            % solve a set of linear equations for each shift
            M = zeros(nTiles, nTiles);
            c = zeros(nTiles, 2);
            for i = 1:size(VerticalPairs,1)
                if isnan(vShifts(i,1)); continue; end
                t1 = VerticalPairs(i,1);
                t2 = VerticalPairs(i,2);
                M(t1,t1) = M(t1,t1)+1;
                M(t1,t2) = M(t1,t2)-1;
                c(t1,:) = c(t1,:) - vShifts(i,:);
                M(t2,t2) = M(t2,t2)+1;
                M(t2,t1) = M(t2,t1)-1;
                c(t2,:) = c(t2,:) + vShifts(i,:);
            end
            
            for i = 1:size(HorizontalPairs,1)
                if isnan(hShifts(i,1)); continue; end
                t1 = HorizontalPairs(i,1);
                t2 = HorizontalPairs(i,2);
                M(t1,t1) = M(t1,t1)+1;
                M(t1,t2) = M(t1,t2)-1;
                c(t1,:) = c(t1,:) - hShifts(i,:);
                M(t2,t2) = M(t2,t2)+1;
                M(t2,t1) = M(t2,t1)-1;
                c(t2,:) = c(t2,:) + hShifts(i,:);
            end
            
            % anchor one of the tiles to a fixed coordinate
            Huge = 1e6;
            TileDistFromCenter = abs(mod(0:nTiles-1, nY)-nY/2) + ...
                abs(floor((0:nTiles-1)/nY)-nX/2);
            [~, HomeTile] = min(TileDistFromCenter(:)./~EmptyTiles(:));
            %sub2ind([nY nX], ceil(nY/2), ceil(nX/2));
            M(nTiles+1,HomeTile) = 1;
            c(nTiles+1,:) = [Huge, Huge];
            
            Tiny = 1e-4; % for regularization
            TileOffset0 = (M+Tiny*eye(nTiles+1, nTiles))\c;
            
            % find tiles that are connected to the home tile
            AlignedOK = (TileOffset0(:,1)>Huge/2);
            TileOffset1 = nan(nTiles, 2);
            TileOffset1(AlignedOK,:) = TileOffset0(AlignedOK,:)-Huge;
            
            % RefPos(t,1:2) is origin of reference tile
            RefPos = bsxfun(@minus, TileOffset1, nanmin(TileOffset1))+1;
            
            % tile origin(t,1:2,r)
            TileOrigin =  round(RefPos);
            csvwrite(fullfile(outdir, ['tile_origin_' output_prefix, '.csv']), TileOrigin);
            
            clear RefImages
        end
    end
    
    % stitch image and write
    for c = 1:nChannels
        MaxTileLoc = max(TileOrigin);
        BigIm = zeros(ceil((MaxTileLoc + 2048)), 'uint16');
        for t = NonemptyTiles
            MyOrigin = TileOrigin(t,:);
            if ~isfinite(MyOrigin(1)); continue; end
            if mod(t,10)==0; fprintf('Loading channel %d tile %d image\n', c, t); end
            LocalIm = imread(fullfile(indir, tilefiles{t}), c);
            BigIm(floor(MyOrigin(1))+(1:2048), floor(MyOrigin(2))+(1:2048)) ...
                = max(cat(3, BigIm(floor(MyOrigin(1))+(1:2048), floor(MyOrigin(2))+(1:2048)), imresize(LocalIm, 1)), [], 3);
        end
        imwrite(BigIm, fullfile(outdir, [output_prefix, '_c', num2str(c) '_stitched.tif']));
    end
    
end
end



%% MIST stitching
function stitch_2d_mist(indir, files, outdir, nChannels, ref_channel, ref_file)

outdir_top = outdir;
if nargin < 6 || ref_file==0
    ref_file = 0;
    outdir = fullfile(outdir, 'Stitched2DTiles_MIST_Individual');
else
    outdir = fullfile(outdir, ['Stitched2DTiles_MIST_Ref' num2str(ref_file)]);
end

mkdir(outdir);

s_order = 1:numel(files);
if ref_file
    s_order = [ref_file, setdiff(s_order, ref_file)];
    
    % transformation from tile alignment
    load(fullfile(outdir_top, 'Aligned2DTiles', 'transformation.mat'));
    TForm = cat(2, TForm{:});
    TForm = reshape(TForm, size(TForm,1), 2, []);
    aligntform = zeros(length(TForm), numel(files), 2);
    aligntform(:,:,1) = squeeze(TForm(:,1,:))';
    aligntform(:,:,2) = squeeze(TForm(:,2,:))';
    % reference image
    aligntform(:,:,1) = bsxfun(@minus, aligntform(:,:,1), aligntform(:,ref_file,1)); 
    aligntform(:,:,2) = bsxfun(@minus, aligntform(:,:,2), aligntform(:,ref_file,2)); 
    
end

save_stitched_image = 1;

imwrite(zeros(2048, 'uint16'), [indir '\empty.tif']);

for s = s_order
    
    imfile = files{s};
    [~, output_prefix] = fileparts(imfile);
    
    xypos = csvread(fullfile(outdir_top, '2DTiles', ['tile_coordinates_' output_prefix '.csv']));
    
    [gridorder, tilefiles] = tilepos2grid(xypos,...
        [output_prefix '_t'], '.tif', 0);
    
    output_prefix = [output_prefix '_'];

    if ref_file && s~=ref_file
        % use exisiting tramsformation
        imfile = files{ref_file};
        [~, output_prefix_ref] = fileparts(imfile);
        output_prefix_ref = [output_prefix_ref '_'];
        
        assemble_from_metadata = 1;     % apply the same to other images
        ref_transform = fullfile(outdir, [output_prefix_ref 'metadata-' num2str(ref_channel) '.mat']);
        
        for c = [ref_channel, setdiff(1:nChannels, ref_channel)]
            copyfile(ref_transform, fullfile(outdir, [output_prefix 'metadata-' num2str(ref_channel) '.mat']));
            
            % modify transformation based on tile alignment
            load(ref_transform, 'global_x_img_pos', 'global_y_img_pos');
            [tilenum, idx] = sort(gridorder(:));
            
            global_x_img_pos(idx(tilenum~=0)) = global_x_img_pos(idx(tilenum~=0)) + aligntform(:,s,1);
            global_y_img_pos(idx(tilenum~=0)) = global_y_img_pos(idx(tilenum~=0)) + aligntform(:,s,2);
            save(fullfile(outdir, [output_prefix 'metadata-' num2str(ref_channel) '.mat']),...
                'global_x_img_pos', 'global_y_img_pos', '-append');
            
            stitch_time_slice(indir, tilefiles,...
                outdir, output_prefix,...
                ref_channel, NaN, NaN, 'Overlay', 0,... 
                save_stitched_image, assemble_from_metadata,...
                fullfile(outdir, ['log_' output_prefix num2str(c) '.txt']),...
                10, 10, c);
            
            % crop/pad upper left part because starting position is always
            % set to (1,1) when assembling images for convenience
            try
                I = imread(fullfile(outdir, [output_prefix 'stitched-' num2str(c) '.tif']));
                I = padimg(I, round(min(global_x_img_pos(:)))-1, round(min(global_y_img_pos(:)))-1, 'NW');
                imwrite(I, fullfile(outdir, [output_prefix 'stitched-' num2str(c) '.tif']));
            catch ME
                if contains(ME.identifier, 'fileDoesNotExist')
                    load(fullfile(outdir, [output_prefix 'stitched-' num2str(c) '.tif.mat']));
                    I = padimg(I, round(min(global_x_img_pos(:)))-1, round(min(global_y_img_pos(:)))-1, 'NW');
                    save(fullfile(outdir, [output_prefix 'stitched-' num2str(c) '.tif.mat']), 'I', '-append');
                end
            end
            clear I

        end
        
    else
        for c = [ref_channel, setdiff(1:nChannels, ref_channel)]
            if c == ref_channel
                % always restitch reference round reference channel
                assemble_from_metadata = 0;
            else
                assemble_from_metadata = 1;
            end
            
            stitch_time_slice(indir, tilefiles,...
                outdir, output_prefix,...
                ref_channel, 7, 2, 'Overlay', 0,...
                save_stitched_image, assemble_from_metadata,...
                fullfile(outdir, ['log_' output_prefix num2str(c) '.txt']),...
                10, 10, c);
            
        end
    end
end
end




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



%% inspect random 10 tile images
function inspect_images(indir, files, nTiles, channels)
randomTiles = randperm(nTiles, 9);
Stitched = ls(indir);
Stitched = cellstr(Stitched(3:end,:));

for t = randomTiles
    f = figure;
    hold on;
    for s = 1:numel(files)
        [~, output_prefix] = fileparts(files{s});
        folders = Stitched(contains(Stitched, output_prefix));
        for c = channels
            I = imread(fullfile(indir, folders{c}, ['tile' num2str(t) '.tif']));
            I = imtophat(I, strel('disk', 3));
            [y,x] = find(I==imdilate(I, strel('disk', 3)) & I>1000);
            plot(x,y,'x'); drawnow;
        end
    end
    title(num2str(t));
    legend
    waitfor(f);
end
end



%% analyze and call bases
function Results = test_analysis(indir, files, outdir, nTiles, channels, taglist)
Stitched = ls(indir);
Stitched = cellstr(Stitched(3:end,:));

code = importdata(taglist, ',');
code = cellfun(@(v) strsplit(v, ','), code, 'UniformOutput', 0);
code = cat(1, code{:});
genes = code(:,end);
code = code(:,1);

sym = symlist;
if length(sym) < length(genes)
    sym = repmat(sym, ceil(length(genes)/length(sym)), 1);
    sym = sym(1:length(genes));
end
sym = [sym; {'ch'}];

Results = cell(nTiles,2);
parfor t = 1:nTiles
    t
    
    basecall = [];
    for s = 1:numel(files)
        [~, output_prefix] = fileparts(files{s});
        folders = Stitched(contains(Stitched, output_prefix));
        
        Istack = [];
        for c = channels
            I = imread(fullfile(indir, folders{c}, ['tile' num2str(t) '.tif']));
            I = imtophat(I, strel('disk', 3));
            Istack = cat(3, Istack, I);
        end
        
        [I, base] = max(Istack, [], 3);
        if s == 1
            idx = find(I==imdilate(I, strel('disk', 2)) & I>100);
            [y,x] = ind2sub(size(I), idx);
        end
        basecall = [basecall, base(idx)];
        
    end
    
    if ~isempty(basecall)
        bases = num2barcode(barcode_mat2num(basecall));
        gene = repmat({'NNNN'}, length(bases), 1);
        try
            gene(ismember(bases, code)) = genes(cellfun(@(v) find(strcmp(v, code)), bases(ismember(bases, code))));
        end
        Results(t,:) = [{[x,y]}, {[bases, gene]}];
    
%         clf; plotall(gene, [x,y], ''); update_legend(gca, [genes; {'NNNN'}], sym); title(num2str(t)); drawnow;
    end
end

tilepos = csvread(fullfile(outdir, 'ReslicedTiles', 'tilepos.csv'));
name = cat(1, Results{:,2});
name = name(:,2);
pos = [];
for t = 1:nTiles
    if ~isempty(Results{t,1})
        pos = [pos; bsxfun(@plus, tilepos(t,:), Results{t,1})];
    end
end
save(fullfile(outdir, 'Decoded.mat'), 'Results', 'name', 'pos');

figure; plotall(name, pos, ''); update_legend(gca, [genes; {'NNNN'}], sym);
    
end





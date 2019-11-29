%% iss_suite stitching, modified from @iss/register
% results less optimal, not yet tested extensively
function stitch_2d_iss(indir, files, outdir, nChannels, ref_channel, corr_threshold, ref_file)
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
tile_stitch_channel = 3;

files= image_files;
outdir=output_directory;
 ref_channel=tile_stitch_channel;
 ref_file=reference_round;

indir=fullfile(output_directory, '2DTiles');

 corr_threshold = [.3 .6];


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
    s=1
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

% ISS image analysis workshop, 170614
% Xiaoyan
% tested on MATLAB R2016b

close all; clear;
cd ICP

%% parameters
% o = iss_options;
% o.InputDirectory = '';
% o.OutputDirectory = '';
% o.CodeFile = 'codebook_unique.csv';
% o.bpLabels = {'A', 'C', 'G', 'T'};
% %
% o.CorrThresh=.3;    
% o.ExtraCodes={};
% o.nExtraRounds=0;
% o.DetectionThresh=600;
% o.IsolationThresh=o.DetectionThresh/5;
% o.IsolationRadius1= 2;
% o.IsolationRadius2=7;
% o.Graphics = 1;
% o.PointCloud=1;
% o.SmoothSize = 1;
% 
% 
% load tilefiles
% o.TileFiles = TileFiles_20x([1:5,8],:,:);
% o.TileFiles = cellfun(@(v) strrep(v, 'G:\In process\170315_161220KI_4-3\Preprocessing\', 'FS_tophat\'), o.TileFiles, 'uni', 0);
% save ..\WS170614 o -append

load ..\WS170614 o

%% preprocessing
% make top-hat filtered tile files (for now a dummy function since this has
% already been done)
% o = iss_extract_and_filter_XQ(o);

%% registration and stitching
rr = o.ReferenceRound;
EmptyTiles = strcmp('', squeeze(o.TileFiles(rr,:,:)));
Tiles = find(~EmptyTiles)';
Tiles = Tiles(:)';

[nY, nX] = size(EmptyTiles);
nTiles = nY*nX;
RefTiles = zeros(o.TileSz, o.TileSz, nY, nX, 'uint16');

for t=Tiles(:)'
    [y,x] = ind2sub([nY nX], t);
    if mod(t,10)==0; fprintf('Loading tile %d\n', t); end
    RefTiles(:,:,t) = imread(o.TileFiles{rr,y,x}, o.AnchorChannel);
%     oNoGraphics = o; oNoGraphics.Graphics=0;
%     RefSpots{t} = iss_detect_spots(RefTiles(:,:,t), oNoGraphics);
%     RefSpots{t} = iss_detect_spots(RefTiles(:,:,t), o);
end

VerticalPairs = ~EmptyTiles(1:end-1,:) & ~EmptyTiles(2:end,:);
HorizontalPairs = ~EmptyTiles(:,1:end-1) & ~EmptyTiles(:,2:end);
nVerticalPairs = sum(VerticalPairs(:));
nHorizontalPairs = sum(HorizontalPairs(:));

% two pairs 
figure;
subplot(211); imshow(RefTiles(:,:,2,1)*50); drawnow;
subplot(212); imshow(RefTiles(:,:,3,1)*50); drawnow;

% middle left and bottom left
i = 2;
[y,x] = ind2sub(size(VerticalPairs), i);
ImRegFft(RefTiles(:,:,y,x), RefTiles(:,:,y+1,x), 's', o.CorrThresh, o.MinSize, 1);
pause()

% do for all
% o = iss_register(o);
load o1 

% stitched anchor image
MaxTileLoc = max(o.RefPos);
BigAnchorIm = zeros(ceil((MaxTileLoc + o.TileSz)), 'uint16');
for t=Tiles
    if ~isfinite(o.RefPos(t,1)); continue; end
    LocalAnchorIm = imread(o.TileFiles{o.ReferenceRound,t}, o.AnchorChannel);
    BigAnchorIm(floor(o.RefPos(t,1))+(1:o.TileSz), ...
        floor(o.RefPos(t,2))+(1:o.TileSz)) ...
        = LocalAnchorIm;
end

% neighboring stitching lines
figure, imshow(BigAnchorIm*50); hold on
for t = Tiles([2,4,5,6,8])
    if ~isfinite(o.RefPos(t,1)); continue; end
    plot([o.RefPos(t,2),o.RefPos(t,2); o.RefPos(t,2), o.RefPos(t,2)+o.TileSz; o.RefPos(t,2)+o.TileSz,o.RefPos(t,2)+o.TileSz;o.RefPos(t,2)+o.TileSz,o.RefPos(t,2)],...
        [o.RefPos(t,1),o.RefPos(t,1)+o.TileSz;  o.RefPos(t,1)+o.TileSz,o.RefPos(t,1)+o.TileSz;o.RefPos(t,1)+o.TileSz,o.RefPos(t,1);o.RefPos(t,1),o.RefPos(t,1)],...
        'y:', 'linewidth', 1.5)
end

%% segmentation
load o2

% basic variables
rr = o.ReferenceRound;
EmptyTiles = strcmp('', squeeze(o.TileFiles(rr,:,:)));
Tiles = find(~EmptyTiles)';

[nY, nX] = size(EmptyTiles);
nTiles = nY*nX;

% middle tile
t = 5;
[y,x] = ind2sub([nY nX], t);
AnchorIm = imread(o.TileFiles{rr,y,x}, o.AnchorChannel);
[RawLocalYX, RawIsolated] = iss_detect_spots(AnchorIm, o);

figure; imshow(AnchorIm*20);
hold on;
plot(RawLocalYX(:,2), RawLocalYX(:,1), 'r+',...
    'markersize', 5, 'linewidth', 1);
plot(RawLocalYX(RawIsolated,2), RawLocalYX(RawIsolated,1), 'yo',...
    'markersize', 5, 'linewidth', 1);

%% from all tiles
% spots in anchor channel on ref round
RawLocalYX = cell(nTiles,1);  % cell array, giving spots in local coordinates
RawIsolated = cell(nTiles,1);
for t=Tiles
    if mod(t,10)==0; fprintf('Detect spots in tile %d\n', t); end;
    [y,x] = ind2sub([nY nX], t);
    AnchorIm = imread(o.TileFiles{rr,y,x}, o.AnchorChannel);
    [RawLocalYX{t}, RawIsolated{t}] = iss_detect_spots(AnchorIm, o);
end
    
% global coordinates
AllIsolated = logical(vertcat(RawIsolated{:}));
nAll = length(AllIsolated);

AllGlobalYX = zeros(nAll,2);
AllLocalYX = zeros(nAll,2);

ind = 1;
for t=Tiles
    MySpots = RawLocalYX{t};
    nMySpots = size(MySpots, 1);
    AllGlobalYX(ind:ind+nMySpots-1,:) = bsxfun(@plus, MySpots, o.RefPos(t,:));
    AllLocalYX(ind:ind+nMySpots-1,:) = MySpots;
    OriginalTile(ind:ind+nMySpots-1) = t;
    ind = ind+nMySpots;
end
figure;
plot(AllGlobalYX(:,2), AllGlobalYX(:,1), '.', 'markersize', 5);
title('All global coords including duplicates');
set(gca, 'YDir', 'reverse');
pause()

% remove duplicates by keeping only spots detected on their home tile
[AllLocalTile, ~] = iss_which_tile(AllGlobalYX, o.RefPos, o.TileSz);
NotDuplicate = (AllLocalTile==OriginalTile');
ndGlobalYX = AllGlobalYX(NotDuplicate,:);
ndLocalYX = AllLocalYX(NotDuplicate,:);
ndIsolated = AllIsolated(NotDuplicate,:);
ndLocalTile = AllLocalTile(NotDuplicate,:);

figure; hold on;
plot(ndGlobalYX(:,2), ndGlobalYX(:,1), '.', 'markersize', 5);
plot(AllGlobalYX(~NotDuplicate,2), AllGlobalYX(~NotDuplicate,1), 'x', 'markersize', 5);
set(gca, 'YDir', 'reverse');
title('Not duplicated and duplicated');
pause()

nnd = sum(NotDuplicate);

% find coordinates of each spot in appropriate tile
ndRoundYX = nan(nnd,2,o.nRounds);
ndRoundTile = nan(nnd,o.nRounds);

for r=1:o.nRounds
    fprintf('Finding coordinates for round %d\n', r);
    % for each spot, find which tile it is in for this round
    
    % this array contains the offset a single tile in round r to ref round
    % the last thing is a way of diagonalizing
    SameTileRelativePos = reshape(o.RelativePos(r,:,1:nTiles+1:nTiles^2), [2, nTiles])';
    
    [ndRoundTile(:,r), ~] = iss_which_tile(ndGlobalYX, o.RefPos-SameTileRelativePos, o.TileSz);
    
    % now shift it. sub2ind avoids making a 3d matrix. Maybe use IndexArrayNan?
    %EachSpotShift = squeeze(o.RelativePos(r,:,sub2ind([nTiles nTiles], ndRoundTile(:,r), ndLocalTile)))';
    if 0; %any(~cellfun(@isempty,o.PcTransform(r,:)))
        for t2=Tiles % tile you are coming from (home on ref round)
            for t1=Tiles % tile you are going to (on other round)
                MySpots = (ndRoundTile(:,r)==t1 & ndLocalTile==t2);
                if ~any(MySpots); continue; end;
                X1 = [ndLocalYX(MySpots,:), ones(sum(MySpots),1)];
                Y = clip(round(X1*o.PcTransform{r,t1,t2}),0,o.TileSz-1);
                ndRoundYX(MySpots,:,r) = Y;
            end
        end
    else
        IndexArray = zeros(4, nnd, 2);
        IndexArray(:,:,1) = [repmat([r 1], nnd, 1), ndRoundTile(:,r), ndLocalTile]';
        IndexArray(:,:,2) = [repmat([r 2], nnd, 1), ndRoundTile(:,r), ndLocalTile]';
        EachSpotShift = IndexArrayNan(o.RelativePos, IndexArray);
        ndRoundYX(:,:,r) = ndLocalYX + EachSpotShift;
    end
    
end

%% point cloud
r = 1; b = 1; % hyb1 base A
t = 5;
MySpots = (ndRoundTile(:,r)==t);
FileName = o.TileFiles{r,t};
TifObj = Tiff(FileName);
TifObj.setDirectory(o.AnchorChannel + b);
BaseIm = TifObj.read();

BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
oNoGraphics = o; oNoGraphics.Graphics=0;
BaseYX = iss_detect_spots(BaseIm,oNoGraphics);

% now loop over all potential home tiles for this one
MyRefTiles = unique(ndLocalTile(ndRoundTile(:,r)==t));
t2 = 5;

MyBaseSpots = (ndRoundTile(:,r)==t & ndLocalTile==t2);
MyLocalYX = ndLocalYX(MyBaseSpots,:);
MyShift0 = o.RelativePos(r,:,t,t2);

[M, error] = PointCloudRegister(BaseYX, MyLocalYX, MyShift0, 3);
MyCorrectedYX = round([MyLocalYX, ones(sum(MyBaseSpots),1)]*M);

% RefIm = RefTiles(:,:,5);
% clf; imshow(RefIm*50, [])
% hold on;
% plot(MyLocalYX(:,2), MyLocalYX(:,1), 'r+');

figure; 
% ax1 = subplot(121);
imshow(BaseIm*70, [])
hold on;
% plot([MyLocalYX(:,2), MyCorrectedYX(:,2)]',...
%     [MyLocalYX(:,1), MyCorrectedYX(:,1)]', '-',...
%     'linewidth', 1, 'color', [.7 .7 .7])
% plot(BaseYX(:,2), BaseYX(:,1), '.');
plot(MyLocalYX(:,2), MyLocalYX(:,1), 'g+', 'linewidth', 1.5, 'markersize', 8);
plot(MyCorrectedYX(:,2), MyCorrectedYX(:,1), 'r+', 'linewidth', 1.5, 'markersize', 8);
% ax2 = subplot(122);
% imshow(BaseIm*70, [])
% hold on;
% plot([MyLocalYX(:,2), MyCorrectedYX(:,2)]',...
%     [MyLocalYX(:,1), MyCorrectedYX(:,1)]', '-',...
%     'linewidth', 1, 'color', [.7 .7 .7])
% linkaxes([ax1 ax2], 'xy');
title('before (green) and after (red) ICP of reference blobs on base2 anchor');
pause()
    
x = meshgrid(min(MyLocalYX(:)):40:max(MyLocalYX(:)));
y = meshgrid(min(MyLocalYX(:)):40:max(MyLocalYX(:)))';
YX2 = [y(:), x(:), ones(numel(y),1)]*M;

figure;
imshow(BaseIm*70, [])
hold on;
% plot(x(:),y(:),'.')
% plot(YX2(:,2),YX2(:,1),'.')
plot([x(:), YX2(:,2)]',...
    [y(:), YX2(:,1)]', '-',...
    'linewidth', 1, 'color', 'c')
axis off
set(gca, 'ydir', 'reverse');

% do for all
% [GlobalYX, SpotColors, Isolated] = iss_find_spots(o);
load iss GlobalYX SpotColors Isolated


%% assign them to cells
% crosstalk cluster assignment
nChans = o.nBP+2;
SpotColors = bsxfun(@rdivide, SpotColors, prctile(SpotColors, o.SpotNorm(1), o.SpotNorm(2)));

BleedMatrix = zeros(o.nBP,o.nBP,o.nRounds); % (Measured, Real, Round)
figure;
set(gcf, 'name', 'Vector clustering', 'units', 'normalized', 'position', [0 0 1 1]); 
for r =1:o.nRounds
    m = squeeze(SpotColors(Isolated,:,r)); % data: nCodes by nBases    
    [Cluster, v, s2] = ScaledKMeans(m, eye(4));
    
    set(gcf, 'name', ['base' num2str(r)]);
    for j = 1:4
        for k = 1:4
            subplot(4,4,(j-1)*4+k);
            gscatter(m(:,j), m(:,k), Cluster);
            xlabel(o.bpLabels{j});
            ylabel(o.bpLabels{k});
            legend off
            drawnow
        end
    end
    pause()
    
    for i=1:4
        BleedMatrix(:,i,r) = v(i,:) * sqrt(s2(i));
    end
end
pause()

% bleedthrough matrix
figure;
for i=1:5
    subplot(2,3,i); 
    imagesc(BleedMatrix(:,:,i)); 
    caxis([0 1]); 
    title(sprintf('Cycle %d', i)); 
    set(gca, 'xtick', 1:4);
    set(gca, 'XTickLabel', {'T', 'G', 'C', 'A'});
    set(gca, 'ytick', 1:4);
    set(gca, 'yTickLabel', {'T', 'G', 'C', 'A'});
    if i==4
        xlabel('Actual')
        ylabel('Measured');
    end
end
subplot(2,3,6);
caxis([0 1]); 
axis off
colormap hot
colorbar
pause()

% now load in the code book and apply bleeds to it
codebook_raw = importdata(o.CodeFile);
CharCode = codebook_raw.textdata(2:end,5);
GeneName = codebook_raw.textdata(2:end,3);
nCodes = size(CharCode,1)-2; % bit of a hack to get rid of Sst and Npy

% create numerical code (e.g. 33244 for CCGAA)
NumericalCode = zeros(nCodes, o.nRounds);
for r = 1:o.nRounds
    NumericalCode(:,r) = codebook_raw.data(1:nCodes,(r-1)*nChans + (3:nChans))*(1:o.nBP)';
end

BledCodes = zeros(nCodes, o.nBP*o.nRounds);
UnbledCodes = zeros(nCodes, o.nBP*o.nRounds);
% make starting point using bleed vectors (means for each base on each day)
for i=1:nCodes
    for r=1:o.nRounds
        BledCodes(i,(1:o.nBP) + (r-1)*o.nBP) = BleedMatrix(:, NumericalCode(i,r), r);
        UnbledCodes(i,NumericalCode(i,r) + (r-1)*o.nBP) = 1;
    end
end
NormBledCodes = bsxfun(@rdivide, BledCodes, sqrt(sum(BledCodes.^2,2)));

figure; 
subplot(121); imagesc(UnbledCodes);
subplot(122); imagesc(NormBledCodes);
pause()

FlatSpotColors = SpotColors(:,:);
Intensity = sqrt(sum(FlatSpotColors.^2,2));
NormFlatSpotColors = bsxfun(@rdivide, FlatSpotColors, Intensity);
SpotScores = NormFlatSpotColors * NormBledCodes';

figure; imagesc(SpotScores);
pause()
set(gca, 'ylim', [.5 10.5]);
xlabel('code');
ylabel('spot');
pause()

[MaxScore, BestCode] = max(SpotScores,[],2);
Codes = CharCode(BestCode);
Genes = GeneName(BestCode);

%% produce output figure
figure;
ShowMe = (MaxScore>.9);
iss_make_figure(o, GlobalYX(ShowMe,:), Genes(ShowMe));

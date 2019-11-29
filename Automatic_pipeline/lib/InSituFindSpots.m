%function [GoodGlobalYX, SpotColors, IsolatedGoodSpots] = InSituFindSpots(o)
% [AllSpotsXY, SpotColors, IsolatedSpots] = InSituFindSpots(Options)
%
% Finds spots in a set of tiles from an in-situ sequencing experiment. This
% should be run after the tiles have been stiched (for the reference round)
% using InSituStitchReferenceTiles and all rounds have been registered
% using InSituRegisterRounds. 
% 
% input structure o contains:
% md: metadata structure from MIST program
% TileToGrid: maps each tile onto the tile grid in this structure
% RoundOffsets: nTimes x 2 x nRounds, comes from round registration
% AnchorChan: anchor channel (default 2)
% ReferenceRound: the one that was aligned (default 2)
% nTiles: # of tiles
% nRounds: # of sequencing rounds (default 5)
% nBP: # of base-pair possibilities (default 4, for life on earth)
% Graphics: 1 to do things in an interactive mode (default 0)
% FilePrefix: the entire path of filename, prior to round #
% FileSuffix: all parts of the filename after tile #
%  + all options in DetectSpots.m 
% ImageScale: used to set threshold, try 9th percentile of image
%
% Filenames are of the form [FilePrefix ROUND _t TILE FileSuffix]
% so if FilePrefix was '/Directory/base', 
% and FileSuffix was '_FS_tophat_stack.tif' then the filename would be
% /Directory/base2_t043_FS_tophat_stack.tif' (for ROUND=2 and TILE=43)
% 
% outputs:
% AllSpotsXY: nSpots by 2 array with coordinates of each spot (even if not
% all colors could be read due to out of frame)
% 
% SpotColors: nSpots x nBases x nRounds array giving unnormalized intensity
% in each channel. nan if it was out of range for a particular round
%
% IsolatedSpots: 1 for spots that around surrounded by emptiness

% if nargin<2; o = struct; end
o = DefaultOpts(o, 'AnchorChan', 2);
o = DefaultOpts(o, 'ReferenceRound', 2);
o = DefaultOpts(o, 'Graphics', 1);
o = DefaultOpts(o, 'nTiles', 186);
o = DefaultOpts(o, 'nRounds', 5);
o = DefaultOpts(o, 'nBP', 4);
o = DefaultOpts(o, 'FilePrefix', 'A:\Dropbox\Dropbox (Neuropixels)\161230_161220KI_3-1\FS_tophat_stack\base');
o = DefaultOpts(o, 'FileSuffix', '_FS_tophat_stack.tif');
o = DefaultOpts(o, 'ImageScale', 600);
o = DefaultOpts(o, 'TileSz', [2048 2048]);

o = DefaultOpts(o, 'DetectionThresh', .5*o.ImageScale);
o = DefaultOpts(o, 'IsolationThresh', .1*o.ImageScale);

% to debug, select this to be less tiles
%DoTiles = 135; % 1:o.nTiles
DoTiles = 1:o.nTiles;

%% variable naming conventions:
% spot subgroups:
% All: Any spot included in any tile (includes duplicates)
% nd: only spots whose anchor coordinate is in its home tile (no duplicates)
% Good: Spots for which could be read for all rounds

% coordinate frames or other info
% LocalYX: relative to home tile on the reference round
% LocalTile: number of home tile on the reference round
% GlobalYX: relative to the stitched image on the reference round
% RoundYX: relative to home tile after registration on each round
% RoundTile: number of home tile after registration on each round
% Isolated: binary number, saying if it is isolated
% SpotColors: the answer:


%% loop through all tiles, finding spots in anchor channel on ref round
RawLocalYX = {}; % cell array, giving spots in local coordinates
RawIsolated = {};
for t=DoTiles
    FileName = sprintf('%s%1d_t%03d%s', o.FilePrefix, o.ReferenceRound, t, o.FileSuffix);
    AnchorIm = imread(FileName, o.AnchorChan);
    [RawLocalYX{t}, RawIsolated{t}] = DetectSpots(AnchorIm, o);
end
    
%% now make array of global coordinates
AllIsolated = vertcat(RawIsolated{:});
nAll = length(AllIsolated);

AllGlobalYX = zeros(nAll,2);
AllLocalYX = zeros(nAll,2);

ind = 1;
for t=DoTiles
    MySpots = RawLocalYX{t};
    nMySpots = size(MySpots, 1);
    AllGlobalYX(ind:ind+nMySpots-1,1) = MySpots(:,1) + o.md.global_y_img_pos(o.TileToGrid(t));
    AllGlobalYX(ind:ind+nMySpots-1,2) = MySpots(:,2) + o.md.global_x_img_pos(o.TileToGrid(t));
    AllLocalYX(ind:ind+nMySpots-1,:) = MySpots;
    OriginalTile(ind:ind+nMySpots-1) = t;
    ind = ind+nMySpots;
end
if o.Graphics
    figure(1001)
    plot(AllGlobalYX(:,2), AllGlobalYX(:,1), '.', 'markersize', 1);
    title('All global coords including duplicates');
    %set(gca, 'YDir', 'reverse');
end

%% now remove duplicates by keeping only spots detected on their home tile

[AllTile, ~] = WhichTileAmIIn(AllGlobalYX, o.ReferenceTileLocations, o.TileSz);
NotDuplicate = (AllTile==OriginalTile');
ndGlobalYX = AllGlobalYX(NotDuplicate,:);
ndLocalYX = AllLocalYX(NotDuplicate,:);
ndIsolated = AllIsolated(NotDuplicate,:);
ndTile = AllTile(NotDuplicate,:);

nnd = sum(NotDuplicate);

if o.Graphics
    figure(1002)
    plot(ndGlobalYX(:,2), ndGlobalYX(:,1), '.', 'markersize', 1);
    title('Global coords without duplicates');
    %set(gca, 'YDir', 'reverse');
end


%% correct global coordinates on each round using the registration for their home tile
GridSz = size(o.md.img_name_grid);
% make array for looking up tile ID from grid position
GridToTile = nan(GridSz);
GridToTile(o.TileToGrid) = 1:o.nTiles;


ndRoundYX = nan(nnd,2,o.nRounds);
ndRoundTile = nan(nnd,o.nRounds);

for t=DoTiles
    
    % find all cells whose anchor spot is in this tile
    MySpots = find(ndTile==t);
    nMySpots = length(MySpots);
    
    % find neighboring tiles
    og = o.TileToGrid(t); % original tile grid index
    [yTile, xTile] = ind2sub(GridSz, og);
    Neighbors = sub2ind(GridSz, yTile, xTile); % start with yourself
    NeighborPos = [0 0]; % positions
    if yTile>1 % add north neighbor
        ng = sub2ind(GridSz, yTile-1, xTile); % neighbor tile grid index
        if isfinite(GridToTile(ng)) % will be NaN if neighbor tile is empty
            Neighbors = [Neighbors; ng]; 
            NeighborPos = [NeighborPos; -o.md.Y1(og), -o.md.X1(og)];
        end
    end
    if xTile>1 % add west neighbor
        ng = sub2ind(GridSz, yTile, xTile-1);
        if isfinite(GridToTile(ng)) % will be NaN if neighbor tile is empty
            Neighbors = [Neighbors; ng];
            NeighborPos = [NeighborPos; -o.md.Y2(og), -o.md.X2(og)];
        end
    end 
    if yTile<GridSz(1) % add south neighbor
        ng = sub2ind(GridSz, yTile+1, xTile);
        if isfinite(GridToTile(ng)) % will be NaN if neighbor tile is empty
            Neighbors = [Neighbors; ng]; 
            NeighborPos = [NeighborPos; o.md.Y1(ng), o.md.X1(ng)];
        end
    end
    if xTile<GridSz(2) % add east neighbor
        ng = sub2ind(GridSz, yTile, xTile+1);
        if isfinite(GridToTile(ng)) % will be NaN if neighbor tile is empty
            Neighbors = [Neighbors; ng]; 
            NeighborPos = [NeighborPos; o.md.Y2(ng), o.md.X2(ng)];
        end
    end
    nNeighbors = length(Neighbors);
    
    % find the appropriate tile and local coord for each round separately
    for r=1:o.nRounds
        % find tile origin for all neighbors in all rounds
        Origin = NeighborPos + o.RoundOffsets(GridToTile(Neighbors), :,r);
        
        % do the computation for this iteration
%        MyCorrected = bsxfun(@plus, ndGlobal(MySpots,:), -o.RoundOffsets(t,:,r));
        [MyNeighbor, MyRoundYX] = WhichTileAmIIn(ndLocalYX(MySpots,:), Origin, o.TileSz);
        HasANeighbor = isfinite(MyNeighbor);
        
        % put it in the big array
        ndRoundYX(MySpots(HasANeighbor),:,r) = MyRoundYX(HasANeighbor,:);
        ndRoundTile(MySpots(HasANeighbor),r) = GridToTile(Neighbors(MyNeighbor(HasANeighbor)));
    end
end

%% get spot colors
%ndSpotColors = nan(nnd, o.nBP, o.nRounds);
for t=1:o.nTiles
    if mod(t,10)==0; fprintf('reading spot colors for tile %d\n', t); end
    for r=1:o.nRounds
        MySpots = (ndRoundTile(:,r)==t);
        if ~any(MySpots); continue; end
        FileName = sprintf('%s%1d_t%03d%s', o.FilePrefix, r, t, o.FileSuffix);
        for b=1:o.nBP
            BaseIm = imread(FileName, o.AnchorChan + b);
            BaseImSm = imfilter(double(BaseIm), fspecial('disk', 1));
            %BaseImSm = imdilate(BaseIm, strel('disk', 2));
            %BaseImSm = BaseIm;
            ndSpotColors(MySpots,b,r) = IndexArrayNan(BaseImSm, 1+ndRoundYX(MySpots,:,r)');
        end
    end
end

%% now find those that were detected in all tiles
Good = all(isfinite(ndSpotColors(:,:)),2);
GoodGlobalYX = ndGlobalYX(Good,:);
GoodSpotColors = ndSpotColors(Good,:,:);
GoodHome = ndTile(Good);
GoodIsolated = ndIsolated(Good);

%% call spots
CodeFile = 'codebook_unique.csv';
NormSpotColors = bsxfun(@rdivide, GoodSpotColors, prctile(GoodSpotColors, 98,1));
[Genes, Codes, MaxScore] = InSituCallSpots(NormSpotColors, find(GoodIsolated), CodeFile);

if o.Graphics
    figure(1)
    ScoreOK = (MaxScore>.85);
    InSituPlot(GoodGlobalYX(ScoreOK,:)/4, Genes(ScoreOK), 'dapi-dsx4.tif');
    %InSituPlot(GoodGlobalYX/4, Genes, 'dapi-dsx4.tif', max(10,(MaxScore-.85)*1000));
    
end

%% sanity check
plsz = 7;
if o.Graphics ==2
    
    GoodRoundYX = ndRoundYX(Good,:,:);
    GoodRoundTile = ndRoundTile(Good,:);
    % which ones to display? Choose those that have different tiles on diff
    % rounds
    % PlotSpots = find(~all(bsxfun(@eq, ndRoundTile(:,1), ndRoundTile),2) & all(isfinite(ndRoundTile),2) & ndIsolated);
    % PlotSpots = 10000:10010;
    
    PlotSpots = find(GoodGlobalYX(:,1)>4045*4 & GoodGlobalYX(:,1)<4050*4 & GoodGlobalYX(:,2)>1585*4 & GoodGlobalYX(:,2)<11587*4);
    
    
    for s=PlotSpots(:)'; %PlotSpots(randperm(length(PlotSpots)))'
        figure(91); clf
        for r=1:o.nRounds
            t=GoodRoundTile(s,r);
            
            y0 = GoodRoundYX(s,1,r);
            x0 = GoodRoundYX(s,2,r);
            y1 = max(1,GoodRoundYX(s,1,r) - plsz);
            y2 = min(o.TileSz(1),GoodRoundYX(s,1,r) + plsz);
            x1 = max(1,GoodRoundYX(s,2,r) - plsz);
            x2 = min(o.TileSz(2),GoodRoundYX(s,2,r) + plsz);
            
            fprintf('Spot %d, round %d, tile %d: y=%d, x=%d\n', s, r, t, GoodRoundYX(s,1,r), GoodRoundYX(s,2,r));

            Ylegs = {'Anchor', 'T', 'G', 'C', 'A'};
            for b=0:o.nBP
                FileName = sprintf('%s%1d_t%03d%s', o.FilePrefix, r, t, o.FileSuffix);
                
                BaseIm = imread(FileName, o.AnchorChan + b, 'PixelRegion', {[y1 y2], [x1 x2]});
                BaseImSm = imfilter(double(BaseIm), fspecial('disk', 3));

                subplot(o.nBP+1, o.nRounds, (b)*o.nRounds + r)
%                imagesc([x1 x2], [y1 y2], BaseImSm); hold on
                imagesc([x1 x2], [y1 y2], BaseIm); hold on
                axis([x0-plsz, x0+plsz, y0-plsz, y0+plsz]);
                plot(xlim, [y0 y0], 'w'); plot([x0 x0], ylim, 'w');
                caxis([0 o.DetectionThresh*2]);
                if r==1; ylabel(Ylegs{b+1}); end
                colorbar;
                
                title(sprintf('Round %d, Base %d, Tile %d', r, b, t));
                drawnow
            end
        end
        figure(92); clf
        imagesc(sq(GoodSpotColors(s,:,:)));
        set(gca, 'ytick', 1:4); set(gca, 'yticklabel', {'T', 'G', 'C', 'A'});
        caxis([0 o.DetectionThresh*2]);
        fprintf('XY = (%d, %d) Called as %s, %s, quality %f\n', ...
            GoodRoundYX(s,2), GoodRoundYX(s,1), Codes{s}, Genes{s}, MaxScore(s));
        pause;
    end
               
                % get subimage

            
end

%% sanity check
if o.Graphics
%     figure(1); 
%     clf; set(gcf, 'color', 'k');
%     set(gca, 'color', 'k');
% %     ImFile = 'A:\Dropbox\Dropbox (Neuropixels)\161230_161220KI_3-1\Stitched\stitchteststitched-2.tif';
%     ImFile = 'dapi-dsx4.tif'
%      StitchedIm = imread(ImFile);
%      imagesc(StitchedIm);
%     hold on

%     figure(1); hold on
%     TileSpots = (GoodHome==135);
    
    %CodeFile = 'A:\Dropbox\Dropbox (Neuropixels)\161230_161220KI_3-1\codebook_unique.csv';


% %     TileGlobal = GoodGlobal(TileSpots,:);
%     uGenes = unique(Genes);
%     for g=uGenes(:)'
%         MySpots = strcmp(Genes, g);
%         plot(GoodGlobalYX(MySpots,2)/4, GoodGlobalYX(MySpots,1/4), '.');
%     end
%     legend(uGenes);
%     ChangeGeneSymbols;
%     axis off
end


    return
%%

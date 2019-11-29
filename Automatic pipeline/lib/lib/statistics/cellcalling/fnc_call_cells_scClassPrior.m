function [CellYX ,pSpotCell ,pCellClass] =...
    fnc_call_cells_scClassPrior(CellMap, gSet, name_inroi, pos_inroi, DapiBoundaries)
% adapted from iss suite by Kenneth Harris
% Xiaoyan, 2018


%% prepare
% % exclude genes in cell typing
% ExcludeGenes = {'Vsnl1'};

% round spot coordianted
SpotYX = round(fliplr(pos_inroi));
SpotGeneName = name_inroi;

%% get info about cells
rp = regionprops(CellMap);
CellYX = fliplr(vertcat(rp.Centroid)); % convert XY to YX
CellArea0 = vertcat(rp.Area); 

MeanCellRadius = mean(sqrt(CellArea0/pi))*.5; % the dapi part is only half of the typical radius
RelCellRadius = [sqrt(CellArea0/pi)/MeanCellRadius; 1]; % but here we want the whole thing

%% single cell expression
% SpotGene(nS): which gene is each spot
% MeanClassExp(nK,nG): mean expression of each gene in each class
% Neighbors(nS, nN): closest neighboring cells for each spot
% D(nS, nN): distance penalty for each of these
% GeneNames(nG): name of each gene
% ClassNames(nK): name of each class

[GeneNames, ~, SpotGeneNo] = unique(SpotGeneName);
TotGeneSpots = accumarray(SpotGeneNo,1);
ClassNames = vertcat(unique(gSet.Class, 'stable'), {'Zero'});

nG = length(GeneNames);
nK = length(ClassNames); % last is zero-expression
nC = size(CellYX,1)+1; % last is misreads
nS = size(SpotYX,1);
nN = 3 + 1; % last is misreads (always a neighbor)

% ClassPrior = [.5*ones(1,nK-1)/nK .5];

ClassDisplayNames = ClassNames;

MeanClassExp = zeros(nK, nG);
gSub = gSet.GeneSubset(GeneNames);
ClassPrior = zeros(1, nK);
for k=1:nK-1 % don't include last since it is zero-expression class
    % get mean expression level, and scale by a factor of 0.2
    % (inefficiency)
    MeanClassExp(k,:) = .2 * mean(gSub.ScaleCell(0).CellSubset(ClassNames{k}).GeneExp,2)';
    
    ClassPrior(k) = gSub.CellSubset(ClassNames{k}).nCells;
end
ClassPrior(nK) = sum(ClassPrior);
ClassPrior = ClassPrior/sum(ClassPrior);

SpotReg = .001; % avoid log 
lMeanClassExp = log(MeanClassExp + SpotReg);   
LogClassPrior = log(ClassPrior);

%% find closest cell
% now find each spot's neighboring cells and distances (nS, nN)
[Neighbors, Dist] = knnsearch(CellYX, SpotYX, 'K', nN);
Neighbors(:,end) = nC; % set last neighbor to misreads

D = -Dist.^2./(2*MeanCellRadius^2) - log(2*pi*MeanCellRadius^2); % don't normalize: bigger cells express more
D(:,end) = log(1e-5); % this is log likelihood of misread

%% reads in a cell get bonus
% any inside cell radius given a bonus
SpotInCell = IndexArrayNan(CellMap, (SpotYX)');
if Neighbors(SpotInCell>0,1)~=SpotInCell(SpotInCell>0)
    error('a spot is in a cell not closest neighbor!');
end
InsideCellBonus = 2;
D(SpotInCell>0, 1) = D(SpotInCell>0, 1) + InsideCellBonus;

%% this is area factor relative to that of the average cell
CellAreaFactor = (exp(-RelCellRadius.^2/2)*(1-exp(InsideCellBonus)) + exp(InsideCellBonus)) ...
    / (exp(-1/2)*(1-exp(InsideCellBonus)) + exp(InsideCellBonus));

%% initialize variables for main loop
pSpotNeighb = zeros(nS, nN); % prob each spot goes to each neighboring cell: last assigned to noise
pCellClass = zeros(nC, nK); % prob each cell goes to each class: last has zero expression

% start a spot in cell it is in, otherwise misread
pSpotNeighb(Neighbors==SpotInCell)=1;
pSpotNeighb(SpotInCell==0,end)=1;

% gammas start off as priors
eSpotGamma = ones(nC, nK, nG);
elSpotGamma = ones(nC, nK, nG)*psi(1); % start with r=1 prior, no spots

eGeneGamma = ones(nG,1); % start with just 1

% this is to check convergence
pSpotNeighbOld = zeros(nS, nN);

%% now main loop
rSpot = 2;  % gamma dist prior for each spot
rGene = 20; % gamma dist shape for scaling whole gene
for iteration=1:100    % max iteration 100
    % CellGeneCount(nC, nG): number of copies of each gene in each cell
    CellGeneCount = zeros(nC,nG);
    for n=1:nN-1
        c = Neighbors(:,n);
        CellGeneCount = CellGeneCount + accumarray([c, SpotGeneNo], pSpotNeighb(:,n), [nC,nG]);
    end
    
    %% call cell gammas

    % eSpotGamma(nC, nK, nG); expected gamma parameter
    % elSpotGamma(nC, nK, nG); expected log gamma parameter
    ScaledMean = CellAreaFactor.*reshape(MeanClassExp,[1 nK nG]);
    eSpotGamma = (rSpot+reshape(CellGeneCount,[nC 1 nG]))./(rSpot + ScaledMean);
    elSpotGamma = psi(rSpot+reshape(CellGeneCount,[nC 1 nG])) - log(rSpot + ScaledMean); % expectation of log gamma

    %% call cells
    % ScaledExp(nC, nK, nG): expected expression under current parameters
    ScaledExp = reshape(MeanClassExp,[1 nK nG]) .* reshape(eGeneGamma,[1 1 nG]) .* CellAreaFactor   + SpotReg;

    % pNegBin(nC, nK, nG): negbin parameter
    pNegBin = ScaledExp ./ (rSpot + ScaledExp);
    
    % wCellClass(nC, nK): summed log likelihoods
    wCellClass = sum(reshape(CellGeneCount,[nC 1 nG]).*log(pNegBin) + rSpot*log(1-pNegBin),3) + LogClassPrior;
    
    % pCellClass(nC, nK): probabilities
    pCellClass = LogLtoP(wCellClass')';

    
    %% call spots
    % wSpotCell(nS, nN)
    aSpotCell = zeros(nS, nN);
    for n=1:nN-1 % don't include misread possibility
        c = Neighbors(:,n);
        aSpotCell(:,n) = sum(pCellClass(c,:) .* lMeanClassExp(:,SpotGeneNo)',2) + ...
            sum(pCellClass(c,:) .* bi(elSpotGamma, c, 1:nK, SpotGeneNo), 2);
    end
    wSpotCell = aSpotCell + D ;
    
    pSpotNeighb = LogLtoP(wSpotCell')';
    MeanProbChanged = max(abs(pSpotNeighb(:)-pSpotNeighbOld(:)));
    fprintf('Iteration %d, mean prob change %f\n', iteration, MeanProbChanged)
    Converged = ( MeanProbChanged<.02); % converges when no probabilities have changed more than this
    pSpotNeighbOld = pSpotNeighb;
    
        %% call gene gammas (etas)
    % to count non-background expression of each gene first compute background
    TotPredictedB = accumarray(SpotGeneNo, pSpotNeighb(:,end), [nG 1]);
    % and total spots in zero cells:
    pCellZero = pCellClass(:,nK); % prob a cell is class zero (nC)
    pSpotZero = sum(pSpotNeighb(:,1:nN-1).*pCellZero(Neighbors(:,1:nN-1)),2); % prob a spot comes from cell class zero (nS)
    TotPredictedZ = accumarray(SpotGeneNo, pSpotZero);
    
    % total counts predicted by all cells of each class (nK, nG)
    ClassTotPredicted = shiftdim(sum(eSpotGamma.*pCellClass.*CellAreaFactor,1),1).*(MeanClassExp + SpotReg);
    % total of each gene (nG): 
    TotPredicted = sum(ClassTotPredicted(1:nK-1,:),1)';
    
    eGeneGamma = (rGene + TotGeneSpots - TotPredictedB - TotPredictedZ)./(rGene + TotPredicted);

    
%     %% diagnostics
%     if ~isempty(o.CellCallShowCenter) && (Converged || o.Graphics==2 || i==o.CellCallMaxIter)
%         figure(3985471)
%         
%         % create a new object to do the plotting in our local coordinate
%         % system
%         o.plot(BackgroundImage,[x0 x1 y0 y1]);
% 
%         % lines to show spots connected to cells
%         [~, BestNeighb] = max(pSpotNeighb,[],2);
%         SpotBestNeighb = bi(Neighbors,(1:nS)',BestNeighb(:));
%         rn = SpotBestNeighb<nC & max(abs(SpotYX-o.CellCallShowCenter),[],2)<o.CellCallShowRad;
%         plot([SpotYX(rn,2) , CellYX(SpotBestNeighb(rn),2)]', ...
%             [SpotYX(rn,1) , CellYX(SpotBestNeighb(rn),1)]', 'Color', [.3 .3 .3]);
%         
%         % text to show best classes
%         [~, BestClass] = max(pCellClass(1:end-1,:),[],2);            
%         text(CellYX(:,2), CellYX(:,1), ClassDisplayNames(BestClass), 'color', 'r', 'fontsize', 6);
%         
%         % zoom in to our axis
%         if ~isempty(o.CellCallShowCenter)
%             axis(reshape(fliplr([o.CellCallShowCenter;o.CellCallShowCenter])...
%                 +[-o.CellCallShowRad; o.CellCallShowRad], 1, 4));
%         end
% 
%         %% diagnostics on an individual cell
%         if ~isempty(o.ExampleCellCenter)
%             
%             [~, MyCell] = min(sum((CellYX-o.ExampleCellCenter).^2,2));
%             fprintf('------------------ Cell %d at %.0f,%.0f: -----------------\n', ...
%                 MyCell, CellYX(MyCell,1), CellYX(MyCell,2));
%             for ss=find((Neighbors(:,1)==MyCell))'
%                 fprintf('Spot %d: %s, with prob %f\n', ss, GeneNames{SpotGeneNo(ss)}, pSpotNeighb(ss,1));
%             end
%             fprintf('-- Total Gene Count --\n');
%             for gg=find(CellGeneCount(MyCell,:)>1e-3) 
%                 fprintf('%s:\t%f\n', GeneNames{gg}, CellGeneCount(MyCell,gg)); 
%             end
%             fprintf('-- Class Posteriors --\n');
%             for cc=find(pCellClass(MyCell,:)>1e-3)
%                 fprintf('%s:\t%e\n', ClassDisplayNames{cc}, pCellClass(MyCell,cc)); 
%             end
%             
%             figure(986543)
%             Myp = squeeze(pNegBin(MyCell,:,:)); % nK by nG
%             WeightMap = CellGeneCount(MyCell,:) .* log(Myp) +  o.rSpot*log(1-Myp);
%             imagesc(WeightMap);
%             set(gca, 'xtick', 1:nG); set(gca, 'XTickLabel', GeneNames); set(gca, 'XTickLabelRotation', 90);
%             set(gca, 'ytick', 1:nK); set(gca, 'yTickLabel', ClassNames);
%             title(sprintf('Cell %d: Contribution of genes to class scores', MyCell));
%             
%             figure(19043765)
%             bar(eGeneGamma);
%             set(gca, 'XTick', 1:nG), set(gca, 'XTickLabel', GeneNames);
%             set(gca, 'XTickLabelRotation', 90);
%             title('Gene efficiencies');
%             grid on
% 
% 
%             [~, TopClasses] = sort(pCellClass(MyCell,:), 'descend');
% %              TopClasses(1) = strmatch('Calb2.Cntnap5a.Rspo3', ClassNames);
% %              TopClasses(2) = strmatch('Cck.Cxcl14.Slc17a8', ClassNames);
%             GeneContrib = WeightMap(TopClasses(1),:) -  WeightMap(TopClasses(2),:);
%             [sorted, order] = sort(GeneContrib);
%             figure (986544);
%             bar(sorted); 
%             set(gca, 'XTick', 1:nG), set(gca, 'XTickLabel', GeneNames(order));
%             set(gca, 'XTickLabelRotation', 90);
%             title(sprintf('Cell %d: Score for class %s vs %s', MyCell, ClassNames{[TopClasses(1), TopClasses(2)]}));
% 
%         end
%         %%
% %          keyboard
%     end
    
    if Converged; break; end
    
end

%% make dense array output
pSpotCell = sparse(repmat(1:nS,1,nN)', Neighbors(:), pSpotNeighb(:));

%% prepare colors for pie charts
% find classes to collapse
CollapseMe = zeros(nK,1);
Colors = zeros(nK,3);
ClassCollapse = ca1_pie_colors(ClassNames);

for i=1:size(ClassCollapse,1)
    ClassList = ClassCollapse{i,1};
    for j=1:length(ClassList)
        MyClasses = strmatch(ClassList{j}, ClassNames);
        if length(MyClasses)==0; continue; end
        CollapseMe(MyClasses)=i;
        Colors(MyClasses,:) = repmat(ClassCollapse{i,3},length(MyClasses),1);
        DisplayName(MyClasses) = ClassCollapse(i,2);
    end
end

nColorWheel = sum(CollapseMe==0);

Colors0 = hsv(ceil(nColorWheel*1.2));
Colors(~CollapseMe,:) = Colors0(1:nColorWheel,:); % last is zero

%% make pie charts
figure(43908765)
% figure
clf; 
set(gcf, 'Color', 'w');
set(gca, 'color', 'w');
hold on

% load(o.CellMapFile, 'RelCellRadius');

PieSize = prctile(sum(pSpotCell(:,1:end-1), 1), 75);
for c=1:nC-1    % exclude last "misread" cell
    pMy = pCellClass(c,:);
    WorthShowing = find(pMy>.1);    % do not show anything below 0.1
    if ~isempty(WorthShowing)

        h = pie(pMy(WorthShowing), repmat({''}, sum(WorthShowing>0)));

        for hi=1:length(h)/2
            hno = (hi*2-1);
    %         Index = find(strcmp(h(i*2).String, NickNames), 1);
            set(h(hno), 'FaceColor', Colors(WorthShowing(hi),:));

            % size based on number of reads
            set(h(hno), 'Xdata', get(h(hno), 'Xdata')*PieSize*sum(pSpotCell(:,c)) + CellYX(c,2));
            set(h(hno), 'Ydata', get(h(hno), 'Ydata')*PieSize*sum(pSpotCell(:,c)) + CellYX(c,1));            
            
%             set(h(hno), 'EdgeAlpha', 0);
            set(h(hno), 'EdgeAlpha', 1, 'LineWidth', .1, 'EdgeColor', [.6 .6 .6]);
        end
    end
    
    if mod(c,2000)==0
        drawnow
    end
end

ClassShown = find(any(pCellClass>0.1,1));
ClassDisplayNameShown = DisplayName(ClassShown);
[uDisplayNames, idx] = unique(ClassDisplayNameShown, 'stable');
nShown = length(uDisplayNames);
for k=1:nShown
    h = text(max(SpotYX(:,2))*1.1 - min(SpotYX(:,2))*.11, min(SpotYX(:,1)) + k*range(SpotYX(:,1))/nShown, DisplayName{ClassShown(idx(k))}, 'fontsize', 8);
    set(h, 'color', Colors(ClassShown(idx(k)),:));
end


%% child functions
    function v = IndexArrayNan(a, i)     
        sza = size(a);
        szi = size(i);
        szv = szi(2:end); % size of output array
        
        % make i into a n by D=d1*d2*... matrix
        i2 = i(:,:);
        
        % get number of dimensions in actual array and index array
        nda = length(sza);
        ndi = szi(1);
        
        % AllInRange is a 1xD matrix saying whether all coordinates are in range
        if ndi==nda
            CoordsInRange = (i2>=1 & bsxfun(@le, i2, sza'));
            AllInRange = all(CoordsInRange,1);
        elseif ndi>nda
            % if they are not equal, check that all trailing dimensions of index array are 1
            CoordsInRange = (i2(1:nda,:)>=1 & bsxfun(@le, i2(1:nda,:), sza'));
            ExtrasOne = (i2(nda+1:end,:)==1);
            AllInRange = all(CoordsInRange,1) & all(ExtrasOne,1);
        else
            error('Index array not enough dimensions');
        end
        
        % create linear index from i2. (God i hate matlab for having to do this cell thing)
        LinIndCell = num2cell(i2(:,AllInRange),2);
        
        % create output
        if length(szv)==1 % now i REALLY hate matlab.
            v = nan(szv,1); % because otherwise it returns a square
        else
            v = nan(szv);
        end
        
        if iscell(a)
            v = num2cell(v);
        end
        
        v(AllInRange) = a(sub2ind(sza, LinIndCell{:}));
    end


    function p = LogLtoP(L)
        % p = LogLtoP(L)
        %
        % given a set of log likelihoods L, robustly compute the posterior
        % probabilities p(i,c) = exp(L(i,c))/ sum_i(exp(L(i,c)));
        % operates column-wise
        
        L1 = bsxfun(@minus,L, max(L,[],1));
        eL = exp(L1);
        p = bsxfun(@rdivide, eL, sum(eL,1));
    end
    
    function Y = bi(X, varargin)
        % A = bi(B, i1, i2, i3, ..., iN)
        %
        % Broadcast indexing.
        %
        % X is an N-dimensional array
        % i1 ... iN are arrays giving the indices for each dimension.
        % If these are all M-dimensional arrays, the output Y is an M-dimensional
        % array, each element of which is given by A at the corresponding members of the
        % index arrays
        %
        % if any of the index arrays have size singleton dimensions, these are
        % expanded.
        %
        % E.g. if X=[1 2; 3 4; 5 6]
        % bi(X, [1 3], [1 2]) = [1 6];
        % bi(X, [1 3], [1;2]) = [1 5; 2 6];
        
        ZeroIndexArray = 0; % make array of full size
        for i=2:nargin
            ZeroIndexArray = ZeroIndexArray .* varargin{i-1}; % auto-broadcasting!!!
        end
        
        inds = cell(nargin-1,1);
        for i=2:nargin
            inds{i-1} = ZeroIndexArray + varargin{i-1}; % auto-broadcasting!!!
        end
        
        Y = X(sub2ind(size(X), inds{:}));
    end

    function c = ca1_pie_colors(OriClassNames)
        % creates a cell array to store in o.ClassCollapse to nicely display the
        % colors of CA1 cells
        
        % how much do individual classes vary compared to the group average
        NoiseSize = .2;
        
        randn('state', 1);
        
        % copied from change_gene_symbols
        pc_or_in =   hsv2rgb([.4 .5 .5]);
        less_active =   hsv2rgb([.3 .2 .7]);
        pc =        hsv2rgb([1/3 1 1]);
        pc2 =       hsv2rgb([.27 1 .7]);
        in_general = hsv2rgb([2/3 1 1]);
        
        sst =   hsv2rgb([.55 1 1]);
        pvalb = hsv2rgb([.7 .8 1]);
        ngf =   hsv2rgb([.85 1 1]);
        cnr1 =  hsv2rgb([ 1 1 1]);
        vip =   hsv2rgb([ .13 1 1]);
        cxcl14= hsv2rgb([.1 1 .6]);
        
        ivy  =   hsv2rgb([.85 .5 .6]);
        
        % dictionary of colors
        d = {{'Sst'}, sst ; {'Pvalb'}, pvalb ; {'Cacna2d1.Lhx6'}, ivy ; ...
            {'Cacna2d1.Ndnf'}, ngf ; {'Ntng1'}, pc_or_in; ...
            {'Cck.Cxcl14'}, cxcl14 ; {'Cck.Lmo1', 'Cck.Calca', 'Cck.Sema5a'}, cnr1; ...
            {'Calb2', 'Vip'}, vip};
        
        % output
        c = cell(0,3);
        
        for i=1:size(d,1)
            MyClasses = [];
            for j=1:length(d{i,1})
                MyClasses = union(MyClasses, strmatch(d{i,1}{j}, OriClassNames));
            end
            for j=1:length(MyClasses)
                cn = OriClassNames(MyClasses(j));
                color = d{i,2}*(1-NoiseSize) + NoiseSize.*rand(1,3);
                color(color<0) = 0; color(color>1) = 1;
                c = vertcat(c, {cn, cn{1}, color});
            end
        end
        
        % manually-defined classes
        c = vertcat(c, {{'PC.CA1'}, 'PC CA1', pc});
        c = vertcat(c, {{'PC.CA2', 'PC.CA3'}, 'PC other', less_active});
        c = vertcat(c, {{'Astro', 'Endo', 'Oligo', 'Eryth', 'Vsmc', 'Microglia', 'Choroid'}...
            , 'Non neuron', [.5 .5 .5]});
        c = vertcat(c, {{'Zero'}, 'Uncalled', [.1 .1 .1]});
    end

end

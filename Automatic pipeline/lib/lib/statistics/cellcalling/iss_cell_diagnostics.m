function iss_cell_diagnostics(MyCell, o, gSet, CellMap, GeneNames, CellGeneCount, eGeneGamma, secondclass)
% iss_cell_diagnostics(MyCell, o, gSet, CellMap, GeneNames, CellGeneCount, eGeneGamma)


%%
ClassNames = vertcat(unique(gSet.Class, 'stable'), {'Zero'});

nG = length(GeneNames);
nK = length(ClassNames); % last is zero-expression

% gene count
fprintf('-- Total Gene Count --\n');
for gg=find(CellGeneCount(MyCell,:)>1e-3)
    fprintf('%s:\t%f\n', GeneNames{gg}, CellGeneCount(MyCell,gg));
end

% class posterior
fprintf('-- Class Posteriors --\n');
if nargin > 7
    for cc=find(o.pCellClass(MyCell,:)>1e-3 | strcmp(secondclass, ClassNames)')
        fprintf('%s:\t%e\n', ClassNames{cc}, o.pCellClass(MyCell,cc));
    end
else
    for cc=find(o.pCellClass(MyCell,:)>1e-3)
        fprintf('%s:\t%e\n', ClassNames{cc}, o.pCellClass(MyCell,cc));
    end
end
    
% negative binomial
rp = regionprops(CellMap);
CellArea0 = vertcat(rp.Area);
MeanCellRadius = mean(sqrt(CellArea0/pi))*.5; % the dapi part is only half of the typical radius
RelCellRadius = [sqrt(CellArea0/pi)/MeanCellRadius; 1]; % but here we want the whole thing
CellAreaFactor = (exp(-RelCellRadius.^2/2)*(1-exp(o.InsideCellBonus)) + exp(o.InsideCellBonus)) ...
    / (exp(-1/2)*(1-exp(o.InsideCellBonus)) + exp(o.InsideCellBonus));

MeanClassExp = zeros(nK, nG);
gSub = gSet.GeneSubset(GeneNames);
ClassPrior = zeros(1, nK);
for k=1:nK-1 % don't include last since it is zero-expression class
    MeanClassExp(k,:) = o.Inefficiency * mean(gSub.ScaleCell(0).CellSubset(ClassNames{k}).GeneExp,2)';
    ClassPrior(k) = gSub.CellSubset(ClassNames{k}).nCells;
end

ScaledExp = reshape(MeanClassExp,[1 nK nG]) .* reshape(eGeneGamma,[1 1 nG]) .* CellAreaFactor   + o.SpotReg;

% pNegBin(nC, nK, nG): negbin parameter
pNegBin = ScaledExp ./ (o.rSpot + ScaledExp);

% heatmap: genes contribution to classes
figure(986543)
Myp = squeeze(pNegBin(MyCell,:,:)); % nK by nG
WeightMap = CellGeneCount(MyCell,:) .* log(Myp) +  o.rSpot*log(1-Myp);
imagesc(WeightMap);
set(gca, 'xtick', 1:nG); set(gca, 'XTickLabel', GeneNames); set(gca, 'XTickLabelRotation', 90);
set(gca, 'ytick', 1:nK); set(gca, 'yTickLabel', ClassNames);
title(sprintf('Cell %d: Contribution of genes to class scores', MyCell));

% barplot: gene efficiencies
figure(19043765)
bar(eGeneGamma);
set(gca, 'XTick', 1:nG), set(gca, 'XTickLabel', GeneNames);
set(gca, 'XTickLabelRotation', 90);
title('Gene efficiencies');
grid on

% barplot: comparison between top two classes
[~, TopClasses] = sort(o.pCellClass(MyCell,:), 'descend');
% TopClasses(1) = strmatch('Calb2.Cntnap5a.Rspo3', ClassNames);
% TopClasses(2) = strmatch('Cacna2d1.Lhx6.Reln', ClassNames);

if nargin > 7
    TopClasses(2) = strmatch(secondclass, ClassNames);
end
GeneContrib = WeightMap(TopClasses(1),:) -  WeightMap(TopClasses(2),:);
[sorted, order] = sort(GeneContrib);
figure (986544);
bar(sorted);
set(gca, 'XTick', 1:nG), set(gca, 'XTickLabel', GeneNames(order));
set(gca, 'XTickLabelRotation', 90);
title(sprintf('Cell %d: Score for class %s vs %s', MyCell, ClassNames{[TopClasses(1), TopClasses(2)]}));


end



%% 
SpotNormPrctile = 98;
nRounds = 5;

%% get intensity
spotColors = importdata('blobs_comb.csv');
tilepos = spotColors.data(mod(spotColors.data(:,1),5)==1,[3, 10:11]);
spotColors = cat(3,...
    spotColors.data(mod(spotColors.data(:,1),5)==1,6:9),...
    spotColors.data(mod(spotColors.data(:,1),5)==2,6:9),...
    spotColors.data(mod(spotColors.data(:,1),5)==3,6:9),...
    spotColors.data(mod(spotColors.data(:,1),5)==4,6:9),...
    spotColors.data(mod(spotColors.data(:,1),5)==0,6:9));

%% global position
startpos = getcsvtilepos('Tiled_wref.csv');
pos = zeros(length(tilepos),2);
for i = 1:length(startpos)
    pos(tilepos(:,1)==i,:) = tilepos(tilepos(:,1)==i,2:3) + startpos(i,2:3);
end

%% what used to correct crosstalk?
% candidateGenes = {'KIAA0652' 'PTEN1' 'OXSM' 'TMEM8A' 'CCDC105'};
taglist = importdata('taglist_mouse120.csv');
taglist = cellfun(@(v) strsplit(v, ','), taglist, 'UniformOutput', 0);
taglist = cat(1, taglist{:});
% correctBarcodes = taglist(contains(taglist(:,2), candidateGenes),1);
% candidateGenes = taglist(contains(taglist(:,2), candidateGenes),2);
codeNum = barcode2num(taglist(:,1)); 
% [candidateGenes, correctBarcodes]
% 
% [~, candidateCodes] = barcode2num(correctBarcodes);
[~, rawCode] = max(spotColors, [], 2);
rawCode = squeeze(rawCode);
rawDecode = rawCode(:,1)*1e4 + rawCode(:,2)*1e3 + rawCode(:,3)*1e2 + rawCode(:,4)*1e1 + rawCode(:,5);
rawDecode = cellfun(@(v) find(v==codeNum), num2cell(rawDecode), 'UniformOutput', 0);
nnnn = cellfun(@isempty, rawDecode);
rawDecode(nnnn) = {0};
rawDecode = cell2mat(rawDecode);
% 
% forCorrection = ismember(rawCode, candidateCodes, 'rows');
% nnz(forCorrection)

% isolated spots
IDX = rangesearch(pos, pos, 15);
forCorrection = cellfun(@(v) length(v)==1, IDX);
nnz(forCorrection)

%% intensity clustering and crosstalk
spotColorsNorm = bsxfun(@rdivide, spotColors, prctile(spotColors, SpotNormPrctile));

BleedMatrix = zeros(4, 4, nRounds); % (Measured, Real, Round)
for r = 1:nRounds
    m = squeeze(spotColorsNorm(forCorrection,:,r)); % data: nCodes by nBases
    
    [Cluster, v, s2] = ScaledKMeans(m, eye(4));
    for i = 1:4
        BleedMatrix(:,i,r) = v(i,:) * sqrt(s2(i));
    end
end

figure;
for i = 1:nRounds
    subplot(2,3,i);
    imagesc(BleedMatrix(:,:,i));
    caxis([0 1]);
    title(sprintf('Round %d', i));
    set(gca, 'xtick', 1:4);
    set(gca, 'XTickLabel', {'AF750' 'AF488' 'Cy3' 'Cy5'});
    set(gca, 'ytick', 1:4);
    set(gca, 'yTickLabel', {'AF750' 'AF488' 'Cy3' 'Cy5'});
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

%% create code matrix
taglist = importdata('taglist_mouse120.csv');
taglist = cellfun(@(v) strsplit(v, ','), taglist, 'UniformOutput', 0);
taglist = cat(1, taglist{:});

[~, code] = barcode2num(taglist(:,1));

idx = bsxfun(@plus, code, 0:4:4*nRounds-1);
idx = [reshape(idx', [], 1),...
    reshape(repmat(1:length(taglist),nRounds,1), [], 1)];

codemat = zeros(length(taglist), 4*nRounds);
idx = sub2ind(size(codemat), idx(:,2), idx(:,1));
codemat(idx) = 1;

%% code matrix with crosstalk
codematBled = zeros(size(codemat, 1), nRounds);

for i = 1:size(codemat,1)
    for r = 1:nRounds
        codematBled(i,4*(r-1)+(1:4)) = BleedMatrix(:,code(i,r),r);
    end
end

codematBledNorm = codematBled ./ sqrt(sum(codematBled.^2,2));

%% barcode matching
spotColorsFlat = spotColorsNorm(:,:);
spotIntensity = sqrt(sum(spotColorsFlat.^2,2));
spotColorsFlatNorm = bsxfun(@rdivide, spotColorsFlat, spotIntensity);
spotScores = spotColorsFlatNorm * codematBledNorm';

[spotScore, bestCode] = max(spotScores, [], 2);

figure,imagesc(spotScores);
set(gca, 'XTick', 1:length(taglist), 'XTickLabel', taglist(:,2), 'XTickLabelRotation', 90);

nnz(bestCode == rawDecode)

%% write
towrite = [taglist(bestCode,:), num2cell(pos), num2cell(spotScore)]'; 
fid = fopen('Decoded_wCorrection.csv', 'w');
fprintf(fid, 'Code,Gene,X,Y,score\n');
fprintf(fid, '%s,%s,%f,%f,%f\n', towrite{:});
fclose(fid);

%% bars
% length(rawDecode)
% nnz(~nnnn)
% nnz(bestCode==rawDecode)

figure, 
subplot(221);
hist(spotScores(~nnnn))
title(['originally expected n=' num2str(nnz(~nnnn))]);
subplot(222);
hist(spotScores(nnnn))
title(['originally NNNN n=' num2str(nnz(nnnn))]);
subplot(223);
hist(spotScores(bestCode==rawDecode))
title(['same reads n=' num2str(nnz(bestCode==rawDecode))]);
subplot(224);
hist(spotScores(bestCode~=rawDecode))
title(['different reads n=' num2str(nnz(bestCode~=rawDecode))]);


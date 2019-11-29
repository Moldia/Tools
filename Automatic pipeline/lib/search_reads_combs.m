function uPositives = search_reads_combs(name, pos, distlim,...
    mincomb, searchname, plotouterbox, col)
% uPositives = search_reads_combs(name, pos, distlim,...
%     mincomb, searchname, plotouterbox, col)
% find where any at least mincomb of different reads occur together
% Xiaoyan, 2017

uNames = unique(name);
if nargin >= 5
    [name, pos] = removereads(name, setdiff(uNames, searchname), pos);
end

% return if no reads left
if isempty(name); uPositives =[]; return; end;

[~, ~, idxName] = unique(name);

% range search
[idxNN, dNN] = rangesearch(pos, pos, distlim);
maxNN = max(cellfun(@length, idxNN));
idxNN = cellfun(@(v) [v, zeros(1, maxNN-length(v))], idxNN, 'uni', 0);
idxNN = reshape([idxNN{:}], maxNN, [])';
dNN = cellfun(@(v) [v, nan(1, maxNN-length(v))], dNN, 'uni', 0);
dNN = reshape([dNN{:}], maxNN, [])';

positives = zeros(numel(name), 6);
d = distlim/10;
if nargin <= 3
    mincomb = 2;
end

% step distance
while d <= distlim
    tempidx = idxNN.*(double(dNN<=d));   
    
    % gene name index of neighbors
    nameNN = tempidx;
    nameNN(nameNN~=0) = idxName(nameNN(nameNN~=0));
    
    % need to have at least mincomb of neighbors
    nNN = sum(logical(nameNN), 2);
    query = find(nNN >= mincomb);
    
    % number of different gene neighbors
    nNNGenes =  cellfun(@(v)...
        length(unique(nameNN(v,:))) - double(ismember(0, nameNN(v,:))),...
        num2cell(query));
    positive = nNNGenes >= mincomb;
    nNNGenes = nNNGenes(positive);
    
    % index of positive queries
    positive = query(positive);    
    
    % remove queries that have already processed and have same number of
    % different gene neighbors
    idx = positives(positive,end)==nNNGenes;
    nNNGenes(idx) = [];
    positive(idx) = [];
    
    % add to potential pool
    for i = 1:length(positive)
        groups = tempidx(positive(i),:);
        groups(groups==0) = [];
        centroid = mean(pos(groups,:), 1);
        positives(positive(i),:) =...
            [centroid, d, pos(positive(i),:), nNNGenes(i)];
    end
    
    d = d + distlim/10;
end

% remove any query without neighbor
positives(sum(positives, 2)==0,:) = [];

if ~isempty(positives)
    % merge close clusters
    idxNNcluster = rangesearch(positives(:,1:2), positives(:,1:2), distlim*1.2);
    alreadyCounted = false(size(positives,1),1);
    uPositives = [];
    
    for i = 1:numel(idxNNcluster)
        idxCluster = idxNNcluster{i};
        
        % skip if the query has been processed
        if alreadyCounted(idxCluster(1)); continue; end
        
        duplicates = positives(idxCluster(~alreadyCounted(idxCluster)),:);
        if size(duplicates,1) > 1
            [~, sortDist] = sort(duplicates(:,3));
            duplicates = duplicates(sortDist,:);
            nReads = duplicates(:,6);
            [~, sortnReads] = sort(nReads, 'descend');
            duplicates = duplicates(sortnReads,:);
            uPositives = [uPositives; mean(duplicates(:,1:2),1), duplicates(1,3:6)];
        else
            uPositives = [uPositives; duplicates(1,:)];
        end
        alreadyCounted(idxCluster) = true;
    end
    
    % each query point takes only the cluster that has shortest distance
    % and most number of differerent reads
    [~, sortDist] = sort(uPositives(:,3));
    uPositives = uPositives(sortDist,:);
    nReads = uPositives(:,6);
    [~, sortnReads] = sort(nReads, 'descend');
    uPositives = uPositives(sortnReads,:);
    [~, idx] = unique(uPositives(:,4:5), 'rows');
    uPositives = uPositives(idx,:);
 
    % visualization
    if nargin <= 5
        plotouterbox = 1;
    end
    
    if nargin <= 6
        col = 'y';
    end
    
    % visualization
    boxsize = distlim*2;
    hold on;
    for i = 1:size(uPositives,1)
        plot([uPositives(i,1)-uPositives(i,3), uPositives(i,1)-uPositives(i,3), uPositives(i,1)+uPositives(i,3), uPositives(i,1)+uPositives(i,3), uPositives(i,1)-uPositives(i,3)],...
            [uPositives(i,2)+uPositives(i,3), uPositives(i,2)-uPositives(i,3), uPositives(i,2)-uPositives(i,3), uPositives(i,2)+uPositives(i,3), uPositives(i,2)+uPositives(i,3)],...
            'color', col, 'linewidth', 1.5);
        if plotouterbox
            plot([uPositives(i,4)-boxsize, uPositives(i,4)-boxsize, uPositives(i,4)+boxsize, uPositives(i,4)+boxsize, uPositives(i,4)-boxsize],...
                [uPositives(i,5)+boxsize, uPositives(i,5)-boxsize, uPositives(i,5)-boxsize, uPositives(i,5)+boxsize, uPositives(i,5)+boxsize],...
                'w', 'linewidth', 2);
        end
    end
else
    uPositives = [];
    disp('No occurence of specified reads combination within given distance.');
end

end

function positives = search_reads_cooccur(name, pos, distlim,...
    searchname, plotouterbox, col)
% positives = search_reads_cooccur(name, pos, distlim,...
%     searchname, plotouterbox, col)
% find where all reads occur together
% Xiaoyan, 2017

uNames = unique(name);
if nargin > 3
    [name, pos] = removereads(name, setdiff(uNames, searchname), pos);
end

% return if no read left
if isempty(name); positives =[]; return; end;

% use the least abundant one as query
[uNames, ~, idxName] = unique(name);
cName = hist(idxName, 1:length(uNames));
[~, idxQ] = sort(cName, 'ascend');
idxQ = idxQ(1);
posquery = pos(idxName==idxQ,:);

% range search
[idxNN, dNN] = rangesearch(pos, posquery, distlim);
maxNN = max(cellfun(@length, idxNN));
idxNN = cellfun(@(v) [v, zeros(1, maxNN-length(v))], idxNN, 'uni', 0);
idxNN = reshape([idxNN{:}], maxNN, [])';
dNN = cellfun(@(v) [v, nan(1, maxNN-length(v))], dNN, 'uni', 0);
dNN = reshape([dNN{:}], maxNN, [])';

% step distance
positives = [];
d = distlim/10;
while d <= distlim
    tempidx = idxNN.*(double(dNN<=d));
    nameNN = tempidx;
    nameNN(nameNN~=0) = idxName(nameNN(nameNN~=0));
    nNN = sum(logical(nameNN),2);
    nNNtemp = find(nNN>=length(uNames));
    nameNN = nameNN(nNNtemp,:);
    positive = false(size(nameNN,1),1);
    for i = 1:size(nameNN,1)
        uniNN = unique(nameNN(i,:));
        positive(i) = length(uniNN(uniNN~=0)) == length(uNames);
    end
    positive = nNNtemp(positive);
    
    for i = 1:length(positive)
        groups = tempidx(positive(i),:);
        groups(groups==0) = [];
        centroid = mean(pos(groups,:),1);
        positives = [positives; centroid, d, posquery(positive(i),:)];
    end
    d = d + distlim/10;
end

% visualization
if ~isempty(positives)
    % merge close clusters
    idxNNcluster = rangesearch(positives(:,4:5), positives(:,4:5), distlim);
    alreadyCounted = false(size(positives,1),1);
    uPositives = [];
    
    for i = 1:numel(idxNNcluster)
        idxCluster = idxNNcluster{i};
        
        % skip if the query has been processed
        if alreadyCounted(idxCluster(1)); continue; end
        
        replicates = positives(idxCluster,:);
        if size(replicates,1) > 1
            [~, sortDist] = sort(replicates(:,3));
            replicates = replicates(sortDist,:);
        end
        uPositives = [uPositives; replicates(1,:)];
        alreadyCounted(idxNNcluster{i}) = true;
    end
    
    [~, unipos] = unique(uPositives(:,4:5), 'rows');
    positives = uPositives(unipos,:);
    
    if nargin <= 4
        plotouterbox = 1;
    end
    
    if nargin <= 5
        col = 'y';
    end
    
    % visualization
    boxsize = distlim*2;
    hold on;
    for i = 1:size(positives,1)
        plot([positives(i,1)-positives(i,3), positives(i,1)-positives(i,3), positives(i,1)+positives(i,3), positives(i,1)+positives(i,3), positives(i,1)-positives(i,3)],...
            [positives(i,2)+positives(i,3), positives(i,2)-positives(i,3), positives(i,2)-positives(i,3), positives(i,2)+positives(i,3), positives(i,2)+positives(i,3)],...
            'color', col);
        if plotouterbox
            plot([positives(i,4)-boxsize, positives(i,4)-boxsize, positives(i,4)+boxsize, positives(i,4)+boxsize, positives(i,4)-boxsize],...
                [positives(i,5)+boxsize, positives(i,5)-boxsize, positives(i,5)-boxsize, positives(i,5)+boxsize, positives(i,5)+boxsize],...
                'w', 'linewidth', 2);
        end
    end
 else
    disp('No co-occurence of all reads within given distance.');
end

end

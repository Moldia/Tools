% bin spatial data into hexagonal bins and output count files
% pool all data from files listed in the batch_input_file
% Xiaoyan, 2017


%% modify here
batch_input_file = 'J:\NextPipeline\files_to_pool.txt';
hexagon_radius = 200;   % in pixel
output_fileprefix = 'pooled_bincounts';

%% do not modify

% which files to pool
fid = fopen(batch_input_file, 'r');
files = textscan(fid, '%s', 'delimiter', '\n');
files = files{:};
fclose(fid);

% start processing
allNames = {};
mapNames = cell(numel(files), 1);
binCounts = cell(numel(files), 1);
binPos = cell(numel(files), 1);
for i = 1:numel(files)
    % load name and coordinates and do hexbin
    [name, pos] = getinsitudata(files{i});
    [name, pos] = removereads(name, 'NNNN', pos);
    [uNames, ~, iName] = unique(name);
    [inbin, bincenters] = hexbin(pos, hexagon_radius);
    
    % count transcripts in each bin
    counts = histcounts2(inbin, iName,...
        [unique(inbin)', max(inbin)+1],...
        1:numel(uNames)+1);
    binCounts{i} = counts;
    binPos{i} = [unique(inbin), bincenters];
    
    % add genes that have not appeared yet
    allNames = [allNames, setdiff(uNames, allNames)'];
        
    % name index in the pool
    iuNames = cellfun(@(v) find(strcmp(v, allNames)), uNames);
    mapNames{i} = iuNames;    
end

% organize into one big matrix
rNames = {};
for i = 1:numel(files)
    rNames = [rNames;...
        catstrnum([strtok(files{i}, '.'), '_r', num2str(hexagon_radius), '_hexbin'], binPos{i}(:,1))];
end

matCounts = zeros(length(rNames), numel(allNames));
nrow = 0;
for i = 1:numel(files)
    matCounts(nrow+1:nrow+size(binPos{i}),mapNames{i}) = binCounts{i};
    nrow = nrow + size(binPos{i});
end
[sortedNames, orderNames] = sort(allNames);
matCounts = matCounts(:,orderNames);

% write count file 
binCountsWrite = [rNames, num2cell(matCounts)]';
fid = fopen([output_fileprefix, '_count.csv'], 'w');
fprintf(fid, 'file_bin,');
fprintf(fid, lineformat('%s', numel(sortedNames)), sortedNames{:});
fprintf(fid, ['%s,' lineformat('%d', numel(sortedNames))], binCountsWrite{:});
fclose(fid);

% write hexbin position file
binPos = cat(1, binPos{:});
binPosWrite = [rNames, num2cell(binPos(:,2:3))]';
fid = fopen([output_fileprefix, '_binpos.csv'], 'w');
fprintf(fid, 'file_bin,bincenter_x,bincenter_y\n');
fprintf(fid, '%s,%d,%d\n', binPosWrite{:});
fclose(fid);


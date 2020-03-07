% bin spatial data into hexagonal bins and output count files
% pool all data from files listed in the batch_input_file
% Xiaoyan, 2017


%% modify here
%batch_input_file = 'G:\DIPG pciseq\files_to_poolORGANOIDS.txt';
batch_input_file = 'G:\Glioblastoma_data\files_to_poolDIPG2.txt';
hexagon_radius =200;   % in pixel
output_fileprefix = 'pooled_bincounts_ORGANOIDS';
distx=50;%90
disty=50;%200
output_directory = 'G:\testtsne';
pcacom= 4;   %Number of pca components. Must be between number of samples and number of genes includeda
hexbin_size = distx*0.35;   % size of plots




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
    [inbin, bincenters,iName] = hexbin_alldistr2(pos, hexagon_radius,distx,disty,iName);
    i
    % count transcripts in each bin
    counts = histcounts2(inbin, iName ,...
        [unique(inbin(~ismember(inbin,0)))', max(inbin)+1],...
        1:numel(uNames)+1);
    binCounts{i} = counts;
    binPos{i} = [unique(inbin(~ismember(inbin,0))), bincenters];
    
    % add genes that have not appeared yet
    allNames = [allNames; setdiff(uNames, allNames)];
        
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

rNames=rNames(1:end); %subsection
matCounts=matCounts(1:end,:); %subsection

% write count file 
binCountsWrite = [rNames, num2cell(matCounts)]';
fid = fopen([output_fileprefix, '_count.csv'], 'w');
fprintf(fid, 'file_bin,');
fprintf(fid, lineformat('%s', numel(sortedNames)), sortedNames{:});
fprintf(fid, ['%s,' lineformat('%d', numel(sortedNames))], binCountsWrite{:});
fclose(fid);


% write hexbin position file
binPos = cat(1, binPos{:}); 
binPos = binPos(1:end,:); %subsection
binPosWrite = [rNames, num2cell(binPos(:,2:3))]';
fid = fopen([output_fileprefix, '_binpos.csv'], 'w');
fprintf(fid, 'file_bin,bincenter_x,bincenter_y\n');
fprintf(fid, '%s,%d,%d\n', binPosWrite{:});
fclose(fid);

mean(sum(matCounts,2))

%cor=corrcoef(matCounts);
%clustergram(cor,'ColumnLabels',uNames,'RowLabels',uNames);

hexbin_counts = [output_fileprefix, '_count.csv'];  % already binned data (pooled or not), requires header, compatible with output from BatchHexBin
hexbin_position = [output_fileprefix, '_binpos.csv'];    % bin position file generated at the same time as bin count file



denseTsne_BatchFiles_complementary_function(hexbin_counts,hexbin_position,hexbin_size,output_directory,pcacom);
% 
% % %CLUSTERING
% eucD = pdist(matCounts,'euclidean');
% clustTreeEuc = linkage(eucD,'average'); 
% cophenet(clustTreeEuc,eucD)
% 
% [h,nodes] = dendrogram(clustTreeEuc,0);
% h_gca = gca;
% h_gca.TickDir = 'out';
% h_gca.TickLength = [.002 0];
% h_gca.XTickLabel = []; 
% 
% 
% 
% hidx = cluster(clustTreeEuc,'criterion','distance','cutoff',500); 
% 
% pointsize=300;
% 
% scatter(binPos(:,3),binPos(:,2),pointsize,hidx,'filled','Marker', 's')
% colormap(hsv)   
% 
% 
% 
% 

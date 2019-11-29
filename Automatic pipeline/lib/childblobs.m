function childblobs(blobfile, parentcell, cellprops, outprefix)
% childblobs(blobfile, parentcell, cellprops, outprefix)
% transform data into single cell profile format
% Xiaoyan, 2017


name = getinsitudata(blobfile);
[uniName, ~, idxName] = unique(name);

[uniCell, ~, idxCell] = unique(parentcell);
matrixCount = hist3([idxCell, idxName],...
    [{1:length(uniCell)}, {1:length(uniName)}]);


% cells with blobs
idxParent = cellfun(@(v) find(cellprops(:,1)==v),...
    num2cell(uniCell), 'uni', 0);
nohit = cellfun(@isempty, idxParent);
idxParent(nohit) = {0};
idxParent = cell2mat(idxParent);

try
    cellswBlobs = [uniCell, cellprops(idxParent,2:end)];    
catch
    nNAN = 1;   % cellID=0 or cellID=-1
    while ~idxParent(nNAN)
        nNAN = nNAN+1;
    end
    cellswBlobs = [uniCell,...
        [nan(nNAN-1,4); cellprops(idxParent(nNAN:end),2:end)]];
end

fid = fopen(['Stitched\Cells_wBlobs_' outprefix '.csv'], 'w');
fprintf(fid,'CellID, metadata_position,area,global_x_pos,global_y_pos\n');
fprintf(fid,'%d,%d,%d,%d,%d\n',reshape(cellswBlobs',[],1));
fclose(fid);

% single cell profile, all included
fid = fopen(['CellBlobs\CellBlobs_' outprefix '_Raw.csv'],'w');
header = [{'CellID'}, uniName', {'centroidX'}, {'centroidY'}];
fmt = lineformat('%s', length(header));
fprintf(fid, fmt, header{:});
mwrite = [uniCell, matrixCount, cellswBlobs(:,end-1:end)];
fmt = lineformat('%d', length(header));
fprintf(fid, fmt, mwrite');
fclose(fid);

% single cell profile, remove NNNN and cellID=0 and cellID=-1
idxNNNN = strcmp(uniName, 'NNNN');
fid = fopen(['CellBlobs\CellBlobs_' outprefix '.csv'],'w');
header = [{'CellID'}, uniName(~idxNNNN)', {'centroidX'}, {'centroidY'}];
fmt = lineformat('%s', length(header));
fprintf(fid, fmt, header{:});
try
    mwrite(1:nNAN-1,:) = [];
end
mwrite(:,find(idxNNNN)+1) = [];
fmt = lineformat('%d', length(header));
fprintf(fid, fmt, mwrite');
fclose(fid);

end


function remove_overlap_blobs(blobfile, centersize, buffersize)
% remove_overlap_blobs(blobfile, centersize, buffersize)
% remove blobs in tile overlap region
% Xiaoyan, 2017


data = importdata(blobfile);
header = data.textdata(1,:);
data = data.data;

col = find(~cellfun(@isempty, strfind(lower(header), 'location')));

incenter = readsinsqr(data(:,col),...
    [buffersize, buffersize, buffersize+centersize, buffersize+centersize]);

data = data(incenter,:);
data(:,col) = data(:,col) - buffersize;

fname = strsplit(blobfile, '.');
fname = fname{1};

fid = fopen([fname '_NoOverlap.csv'], 'w');
fprintf(fid, lineformat('%s', length(header)), header{:});
fprintf(fid, lineformat('%d', length(header)), data');
fclose(fid);
end


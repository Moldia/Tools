function fname = anonymize_gene_calls(infile, text_column_to_anonymize)
% fname = anonymize_gene_calls(infile, text_column_to_anonymize)
% anonymize gene names 
% Xiaoyan, 2018

if nargin < 2
    text_column_to_anonymize = 2;
end

data = importdata(infile);

header = data.textdata(1, :);
colNumeric = cellfun(@isempty, data.textdata(2,:));
name = data.textdata(2:end, text_column_to_anonymize);

% remove NNNNs
isNNNN = strcmp(name, 'NNNN');
name = name(~isNNNN);
data.data = data.data(~isNNNN, :);

% anonymize gene names
[uNames, ~, iName] = unique(name);
uNames_new = catstrnum('gene', 1:numel(uNames))';

header = header([text_column_to_anonymize, find(colNumeric)]);
towrite = [uNames_new(iName), num2cell(data.data)]';

% write to current directory
fname = strfind(infile, filesep);
fname = fullfile(cd, ['anonymized_' infile(fname(end)+1:end)]);

fid = fopen(fname, 'w');
fprintf(fid, lineformat('%s', numel(header)), header{:});
fprintf(fid, ['%s,', lineformat('%d', nnz(colNumeric))], towrite{:});
fclose(fid);

% gene lookup table
towrite = [uNames, uNames_new]';
fname = strfind(infile, filesep);
fname = fullfile(cd, ['anonymized_genelookup_' infile(fname(end)+1:end)]);
fid = fopen(fname, 'w');
fprintf(fid, '%s,%s\n', towrite{:});
fclose(fid);

end


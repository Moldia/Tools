function [countTable, fname] = details_to_count(infile, isSave, whichcolumn)
% [countTable, fname] = details_to_count(infile, isSave, whichcolumn)
% read details.csv file and get simple count data
% Xiaoyan, 2017

if nargin > 2
    names = getinsitudata(infile, whichcolumn);
else
    names = getinsitudata(infile);
end

[uNames, ~, idxName] = unique(names);
cNames = hist(idxName, 1:length(uNames));

% move NNNN to the end
idxNNNN = strcmp(uNames, 'NNNN');
if nnz(idxNNNN)
    cNames = [cNames(~idxNNNN), cNames(idxNNNN)];
    uNames = [uNames(~idxNNNN); {'NNNN'}];
end

countTable = [uNames, num2cell(cNames')]';

% save
if nargin > 1 && isSave
    fname = strsplit(infile, '.csv');
    fname = [fname{1}, '_count.csv'];
    fid = fopen(fname, 'w');
    fprintf(fid, 'Name,Count\n');
    fprintf(fid, '%s,%d\n', countTable{:});
    fclose(fid);
end

function [countTable, fname] = details_to_codecount(infile, isSave)
% [countTable, fname] = details_to_codecount(infile, isSave)
% read details.csv file and get code_n_count data
% Xiaoyan, 2017


names = getinsitudata(infile);
codes = getinsitudata(infile, 1);
[uNames, ~, idxName] = unique(names);
[uCodes, ~, idxCode] = unique(codes);
cCodes = hist(idxCode, 1:length(uCodes));

pairs = unique([idxCode, idxName], 'rows');

countTable = [uCodes, num2cell(cCodes'), cell(length(uCodes),1)]';
for i = 1:length(uCodes)
    idx = pairs(pairs(:,1)==i,2);
    countTable(3,i) = uNames(idx);
end

% save
if nargin > 1 && isSave
    fname = strsplit(infile, '.csv');
    fname = [fname{1}, '_codecount.csv'];
    fid = fopen(fname, 'w');
    fprintf(fid, 'Code,Count,GeneName\n');
    fprintf(fid, '%s,%d,%s\n', countTable{:});
    fclose(fid);
end

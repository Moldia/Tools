function replace_file_columns(inputfile, columns_to_replace, newcontent, outputname)
% replace_file_columns(inputfile, columns_to_replace, newcontent, outputname)
% replace specific columns with 
% Xiaoyan, 2017

data = readtable(inputfile);
header = data.Properties.VariableNames;
data = table2cell(data);

% columns contain characters
charcolumns = false(1,numel(header));
for i = 1:numel(header)
    if ischar(data{1,i})
        charcolumns(i) = 1;
    end
end

% replace columns
for i = 1:numel(columns_to_replace)
    if isnumeric(newcontent{i})
        data(:,columns_to_replace(i)) = num2cell(newcontent{i});
    else
        data(:,columns_to_replace(i)) = newcontent{i};
    end
end

% write output file
data = data';
    
fid = fopen(outputname, 'w');
fprintf(fid, lineformat('%s', numel(header)), header{:});
fmt = repmat({'%d'}, 1, numel(header));
fmt(charcolumns) = {'%s'};
fmt = strcat(fmt, ',');
fmt = [fmt{:}];
fmt = [fmt(1:end-1), '\n'];
fprintf(fid, fmt, data{:});
fclose(fid);

end

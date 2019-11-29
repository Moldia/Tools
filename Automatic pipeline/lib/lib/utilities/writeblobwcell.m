function writeblobwcell(blobfile, parentcell, outprefix)
% writeblobwcell(blobfile, parentcell, outprefix)
% write detailed blob file with parent cell
% Xiaoyan, 2017


blob = importdata(blobfile);
fid = fopen(['ParentCell\', outprefix, '_details_wCell.csv'], 'w');
fmt = lineformat('%s', size(blob.textdata,2)+1);
header = [blob.textdata(1,:), {'Parent_Cell'}];
fprintf(fid, fmt, header{:});
blobwrite = [blob.textdata(2:end,1:(size(blob.textdata,2)-size(blob.data,2))),...
    num2cell([blob.data, double(parentcell)])]';
% fmt = lineformat('%d', length(header)-2);
% fmt = ['%s,%s,', fmt];

fmt = lineformat('%d', size(blob.data,2)+1);
fmt = [repmat('%s,', 1, length(header)-(size(blob.data,2)+1)), fmt];

fprintf(fid, fmt, blobwrite{:});
fclose(fid);


end


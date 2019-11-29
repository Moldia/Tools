function fname = savereadsinbw(blobfile, imbw, imscale, suffix)
% select reads in ROI and save to file
% Xiaoyan, 2017

fid = fopen(blobfile, 'rt');
header = fgets(fid);
header = strrep(header, char(10), '')
% header = strrep(header, newline, '')      % newline only available from R2016b
nCols = numel(strfind(header, ',')) + 1;
fclose(fid);

switch nCols
    case 9      % v3
        code = getinsitudata(blobfile, 1);
        [name, pos, cell, tilepos, general, qscore, align] = ...
            getinsitudata(blobfile, 2, 1, 3, 4, 5, 6, 7);
        name = [code, name];
        blobprop = [pos, cell, tilepos, general, qscore, align];
    case 8      % v2 with parent cell
        code = getinsitudata(blobfile, 1);
        [name, pos, tilepos, general, qscore, cell] = ...
            getinsitudata(blobfile, 2, 1, 3, 4, 5, 6);
        name = [code, name];
        blobprop = [pos, tilepos, general, qscore, cell];
    case 7      % v2 without parent cell
        code = getinsitudata(blobfile, 1);
        [name, pos, tilepos, general, qscore] = ...
            getinsitudata(blobfile, 2, 1, 3, 4, 5);
        name = [code, name];
        blobprop = [pos, tilepos, general, qscore];
    case 3      % most simplified form
        [name, pos] = getinsitudata(blobfile, 1);
        blobprop = pos;
end
        
header = [header, ',inROI'];
    
posScaled = correctcoord(pos, imscale);
inROI = readsinroi(posScaled, imbw);

readsinROI = [name, num2cell([blobprop, inROI])];
readsinROI = readsinROI(inROI~=0,:)';

fname = strsplit(blobfile, '.csv');
if nargin > 3
    fname = [fname{1}, '_', suffix, '.csv'];
else
    fname = [fname{1}, '_ROI.csv'];
end
fid = fopen(fname, 'w');
fprintf(fid, [header, '\n']);
fprintf(fid, [repmat('%s,',1,size(name,2)), lineformat('%d', size(blobprop,2)+1)], readsinROI{:});
fclose(fid);

end

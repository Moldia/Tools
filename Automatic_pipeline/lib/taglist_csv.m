function taglist_csv(tagfile)
% convert .m taglist to .csv
% Xiaoyan, 2017


fid = fopen(tagfile, 'r');
line = fgets(fid);
taglist = {};
symbol = {};
while line ~= -1
    if line(1) == char(39)      % line starts with single quote mark
        line = strsplit(line, char(39));
        taglist = [taglist; line(2)];
        if length(line)>=4
            symbol = [symbol; line(4)];
        end
    end
    line = fgets(fid);
end
fclose(fid);
       

taglist = cellfun(@(v) strsplit(v), taglist, 'uni', 0);
taglist = [cellfun(@(v) v{1}, taglist, 'uni', 0),...
    cellfun(@(v) v{2}, taglist, 'uni', 0)]';

fid = fopen([tagfile(1:end-2) '.csv'], 'w');
if isempty(symbol)
    fprintf(fid, '%s,%s\n', taglist{:});
else
    taglist = [taglist;symbol'];
    fprintf(fid, '%s,%s,%s\n', taglist{:});
end
fclose(fid);

end
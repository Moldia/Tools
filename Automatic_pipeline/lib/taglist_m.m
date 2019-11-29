function taglist_m(tagfile)
% convert .csv taglist to .m
% Xiaoyan, 2017


taglist = importdata(tagfile);
        
taglist = cellfun(@(v) strsplit(v,','), taglist, 'uni', 0);
if length(taglist{1})==2
    taglist = [cellfun(@(v) v{1}, taglist, 'uni', 0),...
        cellfun(@(v) v{2}, taglist, 'uni', 0)];
else
    taglist = [cellfun(@(v) v{1}, taglist, 'uni', 0),...
        cellfun(@(v) v{2}, taglist, 'uni', 0),...
        cellfun(@(v) v{3}, taglist, 'uni', 0)];
end

tags = strcat({char(39)},taglist(:,1),{' '},taglist(:,2),{char(39)});

fid = fopen([strtok(tagfile,'.') '.m'], 'w');
name = strsplit(tagfile,'/');
name = strtok(name(end),'.');
fprintf(fid,'function taglist=%s\n', name{1});
fprintf(fid,'taglist = {\n');
if size(taglist,2)==2
    tags = strcat(tags,{';'});
    fprintf(fid, '%s\n', tags{:});
else
    tags = [tags,...
        strcat({char(39)},taglist(:,3),{char(39)},{';'})]';
    fprintf(fid, '%s,%s\n', tags{:});
end
fprintf(fid,'};');
fclose(fid);

end
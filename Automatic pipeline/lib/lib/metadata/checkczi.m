function ntiles = checkczi(filelist)
% ntiles = checkczi(filelist)
% check if metadata of a .czi file is correct or not
% Xiaoyan, 2017 

bfreader = loci.formats.Memoizer(bfGetReader(), 0);
ntiles = [];
for i = 1:length(filelist)
    bfreader.setId(filelist{i});
    [~, ~, ~, ~, xypos, ~] =  get_ome_tilepos(bfreader);
    ntiles = [ntiles; size(xypos, 1)];    
    bfreader.close()
    sprintf('%s\t %d', filelist{i}, size(xypos, 1))
end

 
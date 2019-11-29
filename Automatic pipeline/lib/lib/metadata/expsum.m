function expsum(expfolder, imdir)
% expsum(expfolder, imdir)
% summarize experimental and imaging information
% Xiaoyan, 2017


section = strsplit(expfolder, '\');
section = section{end};

fid = fopen(['ExpSum_', section, '.csv'], 'w');
fprintf(fid,'Experiment date,\nImaging date ID,');
fprintf(fid,'\nSequencing order,\nSandwich assay targets\n');

% oligo database file
libProbes = importdata('C:\Users\Xiaoyan\OneDrive\worky\allProbes.xlsm');
allProbes = libProbes.textdata.Sheet1;
allProbes = allProbes(2:end,[1,5]);
allPrimers = libProbes.textdata.primers;
allPrimers = allPrimers(2:end,[1,5]);

% number of probes used
write_targets('probesused.csv', allProbes, 'Probe');

% number of primers used
try
    write_targets('primersused.csv', allPrimers, 'Primer');
end


imfiles = cellstr(ls(imdir));
im = cellfun(@(v) strcmp(v(end), 'i'), imfiles);
imfiles = imfiles(im);
im = cellfun(@(v) strfind(v, section), imfiles, 'uni', 0);
im = cellfun(@isempty, im);
imfiles = imfiles(~im);

% image metadata
immeta = cell(5, length(imfiles));
for i = 1:length(imfiles)
    i
    
    reader = loci.formats.Memoizer(bfGetReader(), 0);
    reader.setId([imdir, '\', imfiles{i}]);    
    omeMeta = reader.getMetadataStore();
    reader.close();
    
    % filename
    immeta{1,i} = [imdir, '\', imfiles{i}];
    % microscope
    if strcmp(omeMeta.getDetectorModel(0,0), 'HDCamC11440-22CU')
        immeta{2,i} = 'Herr Nilsson';
    elseif strcmp(omeMeta.getDetectorModel(0,0), 'HDCamC11440-42U')
        immeta{2,i} = 'Pippi';
    else
        immeta{2,i} = 'Unknown';
    end
    % date
    immeta{3,i} = omeMeta.getImageAcquisitionDate(0).getValue();
    % filters
    for j = 1:omeMeta.getFilterSetCount(0)
        immeta{4,i} = [char(immeta{4,i}) char(omeMeta.getChannelName(0,j-1)) ';'];
    end
    immeta{4,i}(end) = [];
    % z stacks
    immeta{5,i} = num2str(omeMeta.getPixelsSizeZ(0).getValue());
end

fprintf(fid, '\n');
rows = {'Image file name', 'Microscope', 'Imaging date', 'Filters', 'Z stacks'};
for i = 1:5
    fprintf(fid, [rows{i}, ',']);
    fprintf(fid, lineformat('%s', length(imfiles)), immeta{i,:});
end
fclose(fid);



function write_targets(fileQuery, libOligos, name)
    [uniName, nProbes] = find_oligo_genes(...
        [expfolder, '\', fileQuery], libOligos);
    fprintf(fid,'\nTotal number of %ss,%d\n', lower(name), sum(nProbes));
    fprintf(fid,'Total number of targets,%d\n', length(uniName));
    fprintf(fid, 'Targets,');
    fprintf(fid, lineformat('%s', length(uniName)), uniName{:});
    fprintf(fid, [name, ' number,']);
    fprintf(fid, lineformat('%d', length(uniName)), nProbes);
end

end





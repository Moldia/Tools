function write_roi_countfile(filename, Coord, uNames, ROI_counts, ROI_area)
% write_roi_countfile(filename, Coord, uNames, ROI_counts, ROI_area)
% write a result file for transcript counting within ROIs
% Xiaoyan, 2018


nROIs = size(Coord,2);
fid = fopen(filename, 'w');

% counts
fprintf(fid, 'GeneName,');
fprintf(fid, lineformat('%s', nROIs), Coord{1,:});

for i = 1:length(uNames)
    fprintf(fid, '%s,', uNames{i});
    fprintf(fid, lineformat('%d', nROIs), ROI_counts(i,:));
end

% ROI area
fprintf(fid, '\n');
fprintf(fid, 'ROI area (in square pixels),');
fprintf(fid, lineformat('%d', nROIs), ROI_area);

fclose(fid);

end


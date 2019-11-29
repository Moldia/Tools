function [ROI_count, ROI_freq, ROI_proportion, ROI_area, In] = ...
    count_in_polygon(nROIs, Coord, pos, uNames, iName, cName)
% count transcripts within ROI polygons
% Xiaoyan, 2015-8-26


ROI_count = [];
ROI_freq = [];
ROI_proportion = [];
ROI_area = [];
In = [];

pos = correctcoord(pos, 1);

for i = 1:nROIs
    poly = Coord(2:3,i);
    poly = ([poly{1};poly{2}])';
    in = inpolygon(pos(:,1), pos(:,2), poly(:,1), poly(:,2));
    
    In = [In,in];
    
    m = hist(iName(in), 1:length(uNames));
    ROI_count = [ROI_count, m'];
    ROI_freq = [ROI_freq, m'/sum(m)];
    ROI_proportion = [ROI_proportion, (m./cName)'];
    ROI_area = [ROI_area, polyarea(poly(:,1), poly(:,2))];
end
In = logical(sum(In,2));    % used in plotting

end

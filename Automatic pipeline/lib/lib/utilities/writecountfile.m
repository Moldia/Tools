function writecountfile(filename, ROI_exist, ROI_name, name_uni,...
    ROI_count, ROI_freq, ROI_proportion, ROI_area, col)
% writecountfile(filename, ROI_exist, ROI_name, name_uni,...
%     ROI_count, ROI_freq, ROI_proportion, ROI_area, col)
% write a result file for transcript counting within ROIs
% Xiaoyan, 2014-12-13


fd = fopen(filename,'w');
ROI_number = length(ROI_exist);

% count
fprintf(fd,'Frequency\n');
fprintf(fd,'%s','GeneName');
for j = 1:ROI_number
    fprintf(fd,',%s',ROI_name{j});
end
fprintf(fd,'\n');

for i = 1:length(name_uni)
    fprintf(fd,'%s',name_uni{i});
    for j = 1:ROI_number
        fprintf(fd,',%d',ROI_count(i,j));
    end
    fprintf(fd,'\n');
end

% frequencies within the ROI
fprintf(fd,'\n');
fprintf(fd,'Relative frequency within ROI\n');
fprintf(fd,'%s','GeneName');
for j = 1:ROI_number
    fprintf(fd,',%s',ROI_name{j});
end
fprintf(fd,'\n');

for i = 1:length(name_uni)
    fprintf(fd,'%s',name_uni{i});
    for j = 1:ROI_number
        fprintf(fd,',%.3f',ROI_freq(i,j));
    end
    fprintf(fd,'\n');
end

% ROI/whole section
fprintf(fd,'\n');
fprintf(fd,'Frequency within ROI/Frequency in the whole section\n');
fprintf(fd,'%s','GeneName');
for j = 1:ROI_number
    fprintf(fd,',%s',ROI_name{j});
end
fprintf(fd,'\n');

for i = 1:length(name_uni)
    fprintf(fd,'%s',name_uni{i});
    for j = 1:ROI_number
        fprintf(fd,',%.3f',ROI_proportion(i,j));
    end
    fprintf(fd,'\n');
end

% ROI area
fprintf(fd,'\n');
fprintf(fd,'ROI area (in square pixels)\n');
fprintf(fd,'%s','');

for j = 1:ROI_number
    fprintf(fd,',%s',ROI_name{j});
end
fprintf(fd,'\n');

for j = 1:ROI_number
    fprintf(fd,',%d',ROI_area(j));
end
fprintf(fd,'\n');

% color specification
fprintf(fd,'\n');
fprintf(fd,'%s','Color');
for j = 1:ROI_number
    fprintf(fd,',%s',col{ROI_exist(j)});
end
fprintf(fd,'\n');

fclose(fd);

end


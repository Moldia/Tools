function parentcell = findparentcell...
    (decodematfile, vartoload, allparentcells)
% parentcell = findparentcell...
%     (decodematfile, vartoload, allparentcells)
% load decoding mat file and get parent cell

filteredblobs = load(decodematfile, vartoload);
filteredblobs = filteredblobs.(vartoload);
parentcell = allparentcells(filteredblobs);

end
